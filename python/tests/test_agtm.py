import warnings

import numpy as np
import pytest

import manydream
from manydream import AGTM, DegenerateModelError, DegenerateModelWarning, ManifoldGraph, MultiManifoldModel, datasets
from manydream.agtm import epsilon_net
from manydream.crawling import crawl


def noisy_circle(rng, R=1.0, n=4000, noise=0.03, d=2):
    X = datasets.circle(n, R, d, rng)
    return X, X + noise * rng.normal(size=X.shape)


@pytest.mark.parametrize("R", [1.0, 2.0, 4.0])
def test_curvature_of_circles_is_inverse_radius(rng, R):
    X, _ = noisy_circle(rng, R, int(4000 * R), 0.0)
    g = crawl(X, 0.2, 1, seed=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = AGTM(g).fit(X + 1e-3 * rng.normal(size=X.shape))
    c = m.edge_curvature()
    assert abs(np.median(c) * R - 1) < 0.05
    assert np.quantile(c, 0.1) * R > 0.7 and np.quantile(c, 0.9) * R < 1.3
    assert abs(np.median(m.edge_curvature(method="second")) * R - 1) < 0.4


def test_epsilon_net_properties(rng):
    X = datasets.torus(20000, random_state=rng)
    g = crawl(X, 0.1, 2, seed=1)
    A = g.adjacency()
    from scipy.sparse import csgraph
    for eps in (1, 2, 3):
        c = epsilon_net(A, eps)
        D = csgraph.dijkstra(A, directed=False, unweighted=True, indices=c)
        assert D.min(0).max() <= eps                       # every node is covered
        Dc = D[:, c] + 100 * np.eye(len(c))
        assert Dc.min() > eps                              # centres are separated


def test_fit_is_monotone_and_density_is_normalised(rng):
    X, T = noisy_circle(rng)
    g = crawl(X, 0.2, 1, seed=0)
    m = AGTM(g, max_iter=40, tol=0).fit(T)
    assert np.all(np.diff(m.loglik_) > -1e-9) and m.loglik_[-1] > m.loglik_[0]
    assert m.diagnostics_["status"] in ("max_iter", "converged")
    assert np.linalg.norm(m.means_ - g.nodes, axis=1).max() < 0.1
    # numerical integral of the density over a grid
    a = np.linspace(-1.6, 1.6, 401)
    G = np.stack(np.meshgrid(a, a), -1).reshape(-1, 2)
    assert abs(np.exp(m.score_samples(G)).sum() * (a[1] - a[0]) ** 2 - 1) < 1e-3
    R = m.responsibilities(T[:200])
    assert R.shape == (200, g.n_nodes) and np.allclose(R.sum(1), 1)
    assert np.allclose(m.node_means(T, np.ones(len(T))), 1)
    nm = m.node_means(T, T)
    assert np.linalg.norm(nm - m.means_, axis=1).max() < 0.05
    S = m.sample(5000, random_state=0)
    assert abs(np.linalg.norm(S, axis=1).mean() - 1) < 0.02
    assert abs(m.score(S) - m.score(T)) < 0.15
    # manifold aligned noise: components are elongated along the circle
    w = np.linalg.eigvalsh(m.covariances_)
    assert np.median(w[:, 0] / w[:, 1]) < 0.5


def test_pruned_estep_matches_full_and_solvers_agree(rng):
    X = datasets.torus(15000, random_state=rng)
    T = X + 0.02 * rng.normal(size=X.shape)
    g = crawl(X, 0.12, 2, seed=0)
    assert g.n_nodes > 64
    fits = {}
    for key, kw in {"pruned": {}, "full": {"prune": None}, "cg": {"solver": "cg"}, "direct": {"solver": "direct"}}.items():
        fits[key] = AGTM(g, max_iter=8, tol=0, **kw).fit(T)
        assert np.all(np.diff(fits[key].loglik_) > -1e-9)
    assert np.allclose(fits["pruned"].loglik_, fits["full"].loglik_, atol=1e-7)
    assert np.allclose(fits["cg"].loglik_, fits["direct"].loglik_, atol=2e-3)
    assert np.allclose(fits["pruned"].score_samples(T[:500]), fits["full"].score_samples(T[:500]), atol=1e-6)
    far = T[:5] + 50.0     # far from every component: falls back to the full sum
    assert np.all(np.isfinite(fits["pruned"].score_samples(far)))


def graph_with_stray_branch(rng):
    X, T = noisy_circle(rng)
    g = crawl(X, 0.2, 1, seed=0)
    k = g.n_nodes
    stray = np.array([[3.0, 3.0], [3.2, 3.0], [3.4, 3.0]])
    edges = np.vstack([g.edges, [[0, k], [k, k + 1], [k + 1, k + 2]]])
    return ManifoldGraph(np.vstack([g.nodes, stray]), edges, 1), T


def test_degenerate_components_are_flagged_and_frozen(rng):
    g, T = graph_with_stray_branch(rng)
    with pytest.warns(DegenerateModelWarning):
        m = AGTM(g).fit(T)
    d = m.diagnostics_
    assert len(d["degenerate"]) >= 1 and d["frozen"].any() and d["events"]
    assert set(d["degenerate"]) <= set(range(g.n_nodes))
    assert np.all(np.isfinite(m.means_)) and np.all(np.isfinite(m.covariances_))
    assert np.all(np.isfinite(m.score_samples(T))) and np.all(np.diff(m.loglik_) > -1e-9)
    with pytest.raises(DegenerateModelError):
        AGTM(g, on_degenerate="raise").fit(T)
    with pytest.warns(DegenerateModelWarning):
        s = AGTM(g, on_degenerate="stop").fit(T)
    assert s.diagnostics_["status"] == "degenerate" and np.all(np.isfinite(s.score_samples(T)))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = AGTM(g, zeta_bounds=(0.5, 2.0)).fit(T)
    assert b.zeta_.min() >= 0.5 and b.zeta_.max() <= 2.0


def test_collapse_is_detected(rng):
    # noise-free data: the transversal variance vanishes and the scales shrink without bound
    X, _ = noisy_circle(rng, noise=0.0)
    g = crawl(X, 0.2, 1, seed=0)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        m = AGTM(g, max_iter=200, tol=0, collapse_ratio=0.05).fit(X)
    assert any(issubclass(x.category, DegenerateModelWarning) for x in w)
    assert m.zeta_.min() >= 0.05 * 0.5 and np.all(np.isfinite(m.loglik_))


def test_input_validation(rng):
    X, T = noisy_circle(rng)
    g = crawl(X, 0.2, 1, seed=0)
    with pytest.raises(ValueError):
        AGTM(g).fit(np.c_[T, T[:, 0]])
    with pytest.raises(RuntimeError):
        AGTM(g).score_samples(T)
    with pytest.raises(ValueError):
        AGTM(g, on_degenerate="ignore")
    with pytest.raises(ValueError):
        AGTM(ManifoldGraph(np.zeros((1, 2)), np.zeros((0, 2), int), 1))


def test_multi_manifold_model(rng):
    Xa, Ta = noisy_circle(rng)
    Xb, Tb = noisy_circle(rng, R=0.5, n=2000)
    Tb = Tb + [3.0, 0.0]
    ma = AGTM(crawl(Xa, 0.2, 1, seed=0)).fit(Ta)
    mb = AGTM(crawl(Xb + [3.0, 0.0], 0.2, 1, seed=0)).fit(Tb)
    mm = MultiManifoldModel([ma, mb])
    assert np.allclose(mm.weights_, [2 / 3, 1 / 3])
    assert (mm.predict(Ta) == 0).mean() > 0.999 and (mm.predict(Tb) == 1).mean() > 0.999
    assert np.allclose(mm.predict_proba(Ta[:50]).sum(1), 1)
    both = np.vstack([Ta, Tb])
    assert np.allclose(mm.score_samples(both),
                       np.logaddexp(np.log(2 / 3) + ma.score_samples(both), np.log(1 / 3) + mb.score_samples(both)))
    bg = MultiManifoldModel([ma, mb], background=0.2, bounds=([-2, -2], [4, 2]))
    assert np.isclose(bg.weights_.sum(), 0.8)
    assert bg.predict(np.array([[1.5, 1.5], [1.0, 0.0]])).tolist() == [-1, 0]
    S = bg.sample(3000, random_state=0)
    assert S.shape == (3000, 2) and 0.1 < (bg.predict(S) == -1).mean() < 0.3
    with pytest.raises(ValueError):
        MultiManifoldModel([ma], background=0.5)


def test_tiny_graphs_and_coarse_bases_do_not_break(rng):
    X = np.c_[np.linspace(0, 1, 400), np.zeros(400)]
    T = X + 0.01 * rng.normal(size=X.shape)
    g = ManifoldGraph(np.array([[0.0, 0], [0.33, 0], [0.66, 0], [1.0, 0]]), [[0, 1], [1, 2], [2, 3]], 1)
    for eps in (1, 2, 5):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m = AGTM(g, epsilon=eps).fit(T)
        assert np.all(np.isfinite(m.covariances_)) and np.all(np.isfinite(m.score_samples(T)))
        assert np.all(np.isfinite(m.edge_curvature()))
