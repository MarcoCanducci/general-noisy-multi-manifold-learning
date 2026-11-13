import numpy as np
import pytest
from scipy.sparse import csgraph

import examm
from examm import datasets
from examm.crawling import crawl, multi_manifold_crawl


def covered(graph, X, r):
    return examm._core.KDTree(X).radius_union(graph.nodes, r).astype(bool).mean()


def test_circle_is_a_closed_chain(rng):
    X = datasets.circle(3000, 1.0, 3, rng)
    g = crawl(X, 0.2, 1, seed=0)
    assert g.n_edges == g.n_nodes and (g.degree == 2).all()
    assert csgraph.connected_components(g.adjacency())[0] == 1
    assert np.allclose(np.linalg.norm(g.nodes, axis=1), 1.0)
    assert 0.05 < g.edge_lengths.min() and g.edge_lengths.max() <= 0.2
    assert g.frames.shape == (g.n_nodes, 1, 3)
    # tangent frames are unit vectors orthogonal to the radius
    assert np.allclose(np.linalg.norm(g.frames[:, 0], axis=1), 1)
    assert np.abs(np.einsum("ij,ij->i", g.frames[:, 0], g.nodes)).max() < 0.15


@pytest.mark.parametrize("name,r", [("torus", 0.1), ("mobius", 0.12)])
def test_surfaces_respect_degree_cap_and_cover_the_manifold(rng, name, r):
    X = datasets.torus(30000, random_state=rng) if name == "torus" else datasets.mobius(30000, rng)
    g = crawl(X, r, 2, seed=5)
    assert g.degree.max() <= 4 and g.degree.min() >= 1
    assert csgraph.connected_components(g.adjacency())[0] == 1
    assert covered(g, X, r) > 0.995
    assert len({tuple(sorted(e)) for e in g.edges.tolist()}) == g.n_edges
    assert (g.edges[:, 0] != g.edges[:, 1]).all()
    assert g.edge_lengths.max() <= r + 1e-12
    # frames are orthonormal and tangent: project the normal of the fitted plane
    F = g.frames
    assert np.allclose(np.einsum("kad,kbd->kab", F, F), np.eye(2), atol=1e-8)
    # no two nodes closer than beta * r
    D = np.linalg.norm(g.nodes[:, None] - g.nodes[None], axis=2) + 10 * np.eye(g.n_nodes)
    assert D.min() > 0.4 * r


def test_three_torus_in_4d_is_one_3d_manifold():
    X = datasets.torus3(grid=26)
    g = crawl(X, 1.1, 3, seed=0)
    assert g.dim == 3 and g.degree.max() <= 6
    assert covered(g, X, 1.1) > 0.99
    assert csgraph.connected_components(g.adjacency())[0] == 1


def test_multi_manifold_crawl_separates_interlocked_toroids():
    ds = datasets.make_three_toroids(scale=0.6, random_state=3)
    res = examm.preprocess.resolve("saf").fit_transform(ds.X, 0.1)
    idx = examm.dimensionality_index(res.diffused, 0.1).index
    sel = idx == 2
    lab = ds.label[res.indices][sel]
    a = multi_manifold_crawl(res.diffused[sel], 0.1, 2, noisy=res.noisy[sel], random_state=0)
    assert len(a) == 3
    seen = set()
    for g in a:
        counts = np.bincount(lab[g.points], minlength=4)[1:]
        assert counts.max() / counts.sum() > 0.99
        seen.add(int(counts.argmax()))
        assert g.degree.max() <= 4
    assert seen == {0, 1, 2}
    b = multi_manifold_crawl(res.diffused[sel], 0.1, 2, noisy=res.noisy[sel], random_state=0)
    assert all(np.array_equal(x.nodes, y.nodes) and np.array_equal(x.edges, y.edges) for x, y in zip(a, b))
    # the faster epsilon-graph alternative finds the same three structures (possibly fragmented)
    e = examm.epsilon_graph(res.diffused[sel], 0.1, 2, noisy=res.noisy[sel], random_state=0)
    big = [g for g in e if g.n_nodes > 100]
    assert len(big) >= 3
    for g in big:
        counts = np.bincount(lab[g.points], minlength=4)[1:]
        assert counts.max() / counts.sum() > 0.99


def test_crawling_terminates_on_noise_and_bad_input(rng):
    X = rng.uniform(size=(3000, 3))
    gs = multi_manifold_crawl(X, 0.02, 1, random_state=0)   # hardly any neighbours: nothing to find
    assert isinstance(gs, list) and len(gs) == 0
    gs = multi_manifold_crawl(X, 0.2, 2, random_state=0)
    assert all(g.degree.max() <= 4 for g in gs)
    with pytest.raises(ValueError):
        crawl(X, 0.1, 4)
    with pytest.raises(ValueError):
        crawl(X, -1.0, 1)
    with pytest.raises(ValueError):
        crawl(X, 0.1, 1, eta=0.3, beta=0.5)


def test_graph_helpers(rng):
    X = datasets.circle(2000, 1.0, 2, rng)
    g = crawl(X, 0.2, 1, seed=0)
    two = examm.ManifoldGraph(np.vstack([g.nodes, g.nodes + 5]), np.vstack([g.edges, g.edges + g.n_nodes]), 1)
    parts = two.components()
    assert len(parts) == 2 and all(p.n_nodes == g.n_nodes and p.n_edges == g.n_edges for p in parts)
    assert np.isclose(g.total_length(), g.edge_lengths.sum())
    assert abs(g.total_length() - 2 * np.pi) < 0.2
