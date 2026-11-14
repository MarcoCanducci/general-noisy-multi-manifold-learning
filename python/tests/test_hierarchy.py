import warnings

import numpy as np
import pytest

import manydream
from manydream import AGTM, HierarchicalAGTM, ManifoldGraph, build_hierarchy, datasets
from manydream.crawling import crawl


def two_circles(rng, n_out=150):
    a = datasets.circle(4000, 1.0, 2, rng)
    b = datasets.circle(2000, 0.5, 2, rng) + [3.0, 0.0]
    clean = np.vstack([a, b])
    noisy = clean + 0.03 * rng.normal(size=clean.shape)
    out = rng.uniform([-2.5, -2.5], [5.0, 2.5], size=(n_out, 2))
    out = out[(np.abs(np.linalg.norm(out, axis=1) - 1) > 0.25) & (np.abs(np.linalg.norm(out - [3, 0], axis=1) - 0.5) > 0.25)]
    X = np.vstack([noisy, out])
    truth = np.r_[np.zeros(len(noisy), bool), np.ones(len(out), bool)]
    graphs = [crawl(a, 0.2, 1, seed=0), crawl(b, 0.2, 1, seed=0)]
    return X, truth, graphs


def test_weighted_fit_equals_fit_on_repeated_points(rng):
    X = datasets.circle(1500, 1.0, 2, rng)
    T = X + 0.03 * rng.normal(size=X.shape)
    g = crawl(X, 0.2, 1, seed=0)
    w = rng.integers(1, 4, len(T)).astype(float)
    a = AGTM(g, max_iter=6, tol=0).fit(T, sample_weight=w)
    b = AGTM(g, max_iter=6, tol=0).fit(np.repeat(T, w.astype(int), axis=0))
    assert np.allclose(a.loglik_, b.loglik_, atol=1e-8) and np.allclose(a.means_, b.means_, atol=1e-8)
    assert np.isclose(a.n_samples_, w.sum())
    with pytest.raises(ValueError):
        AGTM(g).fit(T, sample_weight=-w)


def test_one_level_hierarchy_density_and_outliers(rng):
    X, truth, graphs = two_circles(rng)
    H = HierarchicalAGTM(background=0.02).set_data(X)
    nodes = H.add_children(graphs)
    assert len(nodes) == 2 and [k for _, k, _ in H.leaves()] == ["model", "model", "background"]
    assert np.isclose(sum(w for _, _, w in H.leaves()), 1.0)
    assert np.all(np.diff(nodes[0].history) > -1e-8)               # restricted likelihood is monotone
    # mixing coefficients follow the data: 4000 : 2000 : outliers
    assert abs(nodes[0].pi - 4000 / len(X)) < 0.02 and abs(nodes[1].pi - 2000 / len(X)) < 0.02
    assert abs(H.root_.pi_residual - truth.mean()) < 0.01
    # responsibilities of each node sum to its share of the data
    assert np.isclose(nodes[0].resp.sum() / len(X), nodes[0].pi, atol=1e-6)
    # the density integrates to one
    lo, hi = H.bounds_
    gx, gy = np.linspace(lo[0], hi[0], 500), np.linspace(lo[1], hi[1], 400)
    G = np.stack(np.meshgrid(gx, gy), -1).reshape(-1, 2)
    assert abs(np.exp(H.score_samples(G)).sum() * (gx[1] - gx[0]) * (gy[1] - gy[0]) - 1) < 0.01
    # outliers
    flagged = H.outliers(X, alpha=1e-3)
    assert flagged[truth].mean() > 0.9 and flagged[~truth].mean() < 0.01
    assert (H.background_probability(X)[truth] > 0.5).mean() > 0.9
    assert np.isclose(H.outliers(X, alpha=None, percentile=5).mean(), 0.05, atol=0.002)
    lab = H.predict(X[~truth])
    assert (lab[:4000] == 0).mean() > 0.99 and (lab[4000:] == 1).mean() > 0.99
    S = H.sample(2000, random_state=0)
    assert S.shape == (2000, 2) and "background" in H.summary()
    with pytest.raises(ValueError):
        H.add_children(graphs)


def test_children_refine_a_coarse_parent(rng):
    # two parallel noisy lines 0.4 apart; the parent is one coarse line between them
    x = rng.uniform(0, 4, 3000)
    up = np.c_[x[:1500], np.full(1500, 0.2)]
    down = np.c_[x[1500:], np.full(1500, -0.2)]
    X = np.vstack([up, down]) + 0.02 * rng.normal(size=(3000, 2))
    t = np.linspace(0, 4, 21)
    coarse = ManifoldGraph(np.c_[t, np.zeros(21)], [[i, i + 1] for i in range(20)], 1)
    fine = [ManifoldGraph(np.c_[t, np.full(21, s)], [[i, i + 1] for i in range(20)], 1) for s in (0.2, -0.2)]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        H = HierarchicalAGTM(background=0.0).set_data(X)
        (parent,) = H.add_children([coarse])
        before = H.score(X)
        kids = H.add_children(fine, parent)
    assert [n.level for n in H.nodes()] == [2, 3, 3] and len(H.leaves()) == 3
    assert H.score(X) > before + 0.5                              # the finer level fits much better
    assert np.isclose(kids[0].pi + kids[1].pi + parent.pi_residual, 1.0)
    assert abs(kids[0].pi - 0.5) < 0.03 and parent.pi_residual < 0.02
    # eq. 30: a child's responsibility never exceeds its parent's
    assert np.all(kids[0].resp <= parent.resp + 1e-12)
    lab = H.predict(X)
    assert (lab[:1500] == 0).mean() > 0.99 and (lab[1500:] == 1).mean() > 0.99
    assert H.residual_probability(X).mean() < 0.02
    assert np.isclose(sum(w for _, _, w in H.leaves()), 1.0)
    # without residual leaves the parent vanishes from the density
    H2 = HierarchicalAGTM(background=0.0, residual=False).set_data(X)
    (p2,) = H2.add_children([coarse])
    H2.add_children(fine, p2)
    assert [k for _, k, _ in H2.leaves()] == ["model", "model"]


def test_build_hierarchy_by_scale(rng):
    X, truth, _ = two_circles(rng)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        H = build_hierarchy(X, [0.2], preprocess=None, dims=[1], random_state=0, background=0.02)
    assert len(H.root_.children) == 2 and len(H.detections_) == 1
    flagged = H.outliers(X)
    assert flagged[truth].mean() > 0.85 and flagged[~truth].mean() < 0.02
    with pytest.raises(RuntimeError):
        HierarchicalAGTM().leaves()


def test_regional_orthogonal_outliers(rng):
    # a circle sampled 20 times more densely on one half than on the other, in 3-D
    th = np.r_[rng.uniform(0, np.pi, 8000), rng.uniform(np.pi, 2 * np.pi, 400)]
    clean = np.c_[np.cos(th), np.sin(th), np.zeros(len(th))]
    X = clean + 0.02 * rng.normal(size=clean.shape)
    # planted points: on the sparse AND the dense half, pushed off the curve across it
    tp = np.array([0.5, 1.5, 2.5, 3.6, 4.6, 5.6])
    planted = np.c_[1.15 * np.cos(tp), 1.15 * np.sin(tp), np.full(6, 0.1)]
    X = np.vstack([X, planted])
    is_planted = np.r_[np.zeros(len(clean), bool), np.ones(6, bool)]
    g = crawl(clean, 0.2, 1, seed=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        H = HierarchicalAGTM(background=0.0).set_data(X)
        (node,) = H.add_children([g])
    m = node.model
    F = m.tangent_frames()
    assert F.shape == (g.n_nodes, 1, 3) and np.allclose(np.linalg.norm(F[:, 0], axis=1), 1)
    radial = m.means_ / np.linalg.norm(m.means_, axis=1, keepdims=True)
    assert np.abs(np.einsum("kd,kd->k", F[:, 0], radial)).max() < 0.2       # tangents are tangent
    res = H.regional_deviation(X)
    assert set(res) >= {"leaf", "node", "distance", "d2", "p_value", "regional_rank", "leaf_name"}
    assert (res["dof"] == 2).all() and (res["leaf"] == 0).all()
    # the orthogonal distance is the distance to the fitted curve (its piecewise-linear skeleton) ...
    V, E = m.means_, g.edges
    a, b = V[E[:, 0]], V[E[:, 1]]
    sub = np.r_[np.flatnonzero(is_planted), np.arange(0, len(clean), 20)]
    P = X[sub]
    t = np.clip(np.einsum("ned,ed->ne", P[:, None] - a[None], b - a) / ((b - a) ** 2).sum(1), 0, 1)
    brute = np.linalg.norm(P[:, None] - (a[None] + t[:, :, None] * (b - a)[None]), axis=2).min(1)
    assert np.abs(res["distance"][sub] - brute).max() < 0.01 and np.corrcoef(res["distance"][sub], brute)[0, 1] > 0.99
    # ... which tracks the distance from the true circle, wherever along it the point sits
    true_off = np.hypot(np.linalg.norm(X[:, :2], axis=1) - 1, X[:, 2])
    assert np.corrcoef(res["distance"], true_off)[0, 1] > 0.75
    assert res["distance"][is_planted].min() > 0.12
    # planted points stand out regionally on both halves; ordinary points of the sparse half do not
    assert (res["regional_rank"][is_planted] > 0.99).all() and (res["p_value"][is_planted] < 1e-3).all()
    sparse_half = (~is_planted) & (X[:, 1] < -0.1)
    assert (res["regional_rank"][sparse_half] > 0.99).mean() < 0.03
    assert np.median(res["p_value"][sparse_half]) > 0.1
    flagged = H.regional_outliers(X, rank=0.99, p_value=1e-3)
    assert flagged[is_planted].all() and flagged[~is_planted].mean() < 0.005
    # the density score, by contrast, penalises the whole sparse half
    s = H.score_samples(X)
    assert np.median(s[sparse_half]) < np.median(s[(~is_planted) & (X[:, 1] > 0.1)]) - 0.5
    # positions along the curve are right
    ang_node = np.arctan2(m.means_[res["node"], 1], m.means_[res["node"], 0])
    ang = np.arctan2(X[:, 1], X[:, 0])
    assert np.abs(np.angle(np.exp(1j * (ang - ang_node)))).max() < 0.25
