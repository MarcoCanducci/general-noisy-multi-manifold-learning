import numpy as np
import pytest

from examm import _core


@pytest.mark.parametrize("d", [1, 2, 3, 5, 8])
def test_radius_queries_match_brute_force(rng, d):
    X = rng.uniform(size=(600, d))
    Q = rng.uniform(size=(80, d))
    r = 0.25 * d**0.5
    tree = _core.KDTree(X)
    D = np.linalg.norm(Q[:, None] - X[None], axis=2)
    assert np.array_equal(tree.radius_count(Q, r), (D <= r).sum(1))
    indptr, ind, dist = tree.radius_neighbors(Q, r)
    for i in range(len(Q)):
        got = ind[indptr[i]:indptr[i + 1]]
        assert np.array_equal(got, np.flatnonzero(D[i] <= r))
        assert np.allclose(dist[indptr[i]:indptr[i + 1]], D[i, got])
    flags = (rng.uniform(size=len(X)) < 0.1).astype(np.uint8)
    assert np.array_equal(tree.radius_any(Q, r, flags).astype(bool), ((D <= r) & (flags > 0)).any(1))
    assert np.array_equal(tree.radius_union(Q, r).astype(bool), (D <= r).any(0))


@pytest.mark.parametrize("d", [2, 3, 6])
def test_nearest_with_mask(rng, d):
    X = rng.normal(size=(500, d))
    Q = rng.normal(size=(60, d))
    alive = (rng.uniform(size=len(X)) < 0.5).astype(np.uint8)
    tree = _core.KDTree(X)
    D = np.linalg.norm(Q[:, None] - X[None], axis=2)
    idx, dist = tree.nearest(Q)
    assert np.array_equal(idx, D.argmin(1)) and np.allclose(dist, D.min(1))
    Dm = np.where(alive > 0, D, np.inf)
    idx, dist = tree.nearest(Q, 0.8, alive)
    exp = np.where(Dm.min(1) <= 0.8, Dm.argmin(1), -1)
    assert np.array_equal(idx, exp)


@pytest.mark.parametrize("d", [2, 3, 5, 8])
def test_local_pca_matches_numpy(rng, d):
    X = rng.normal(size=(400, d)) * np.linspace(1, 0.2, d)
    tree = _core.KDTree(X)
    r = 0.8 * d**0.5
    counts, ev, vec = tree.local_pca(X[:50], r, d, 2)
    for i in range(50):
        nb = X[np.linalg.norm(X - X[i], axis=1) <= r]
        assert counts[i] == len(nb)
        w, v = np.linalg.eigh(np.cov(nb.T).reshape(d, d))
        assert np.allclose(ev[i], w[::-1], atol=1e-10)
        C = np.cov(nb.T).reshape(d, d)
        assert np.allclose(vec[i] @ vec[i].T, np.eye(d), atol=1e-9)
        assert np.allclose(vec[i] @ C @ vec[i].T, np.diag(w[::-1]), atol=1e-9)


def test_duplicates_and_degenerate_input():
    X = np.zeros((100, 3))
    tree = _core.KDTree(X)
    assert tree.radius_count(X[:1], 0.1)[0] == 100
    with pytest.raises(ValueError):
        _core.KDTree(np.array([[np.nan, 0.0]]))
    with pytest.raises(ValueError):
        tree.radius_count(np.zeros((3, 2)), 0.1)


def test_fortran_ordered_input_is_handled(rng):
    X = np.asfortranarray(rng.uniform(size=(200, 3)))
    a = _core.KDTree(X).radius_count(X, 0.2)
    b = _core.KDTree(np.ascontiguousarray(X)).radius_count(np.ascontiguousarray(X), 0.2)
    assert np.array_equal(a, b)


def test_weighted_cov_and_estep(rng):
    d, K, N = 3, 7, 300
    T = rng.normal(size=(N, d))
    C = rng.normal(size=(K, d))
    scales = rng.uniform(0.5, 1.0, K)
    tree = _core.KDTree(T)
    cov, mass = tree.weighted_cov(C, scales, 3.0, 50.0)
    for k in range(K):
        diff = T - C[k]
        w = np.exp(-(diff**2).sum(1) / (3 * scales[k] ** 2))
        assert np.allclose(cov[k], (diff.T * w) @ diff / w.sum())
        assert np.isclose(mass[k], w.sum())
    # E-step against a direct evaluation
    from scipy.stats import multivariate_normal
    covs = cov + 0.05 * np.eye(d)
    L = np.linalg.cholesky(covs)
    wh = np.linalg.inv(L)
    logdet = 2 * np.log(np.diagonal(L, axis1=1, axis2=2)).sum(1)
    logc = np.log(1.0 / K) - 0.5 * logdet - 0.5 * d * np.log(2 * np.pi)
    dens = np.stack([multivariate_normal(C[k], covs[k]).pdf(T) for k in range(K)], 1) / K
    R = dens / dens.sum(1, keepdims=True)
    vals = rng.normal(size=(N, 2))
    logp, G, S, CC, Vs = _core.agtm_estep(T, C, wh, logc, vals, True)
    assert np.allclose(logp, np.log(dens.sum(1)))
    assert np.allclose(G, R.sum(0)) and np.allclose(S, R.T @ T) and np.allclose(Vs, R.T @ vals)
    assert np.allclose(CC, np.einsum("nk,na,nb->kab", R, T, T))
    assert np.allclose(_core.agtm_resp(T, C, wh, logc), R)
