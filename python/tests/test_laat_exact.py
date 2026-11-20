import numpy as np
import pytest

from manydream import laat_exact
from manydream.laat_exact import LAATGraph, run, walk_lengths


@pytest.fixture(scope="module")
def data():
    g = np.random.default_rng(0)
    t = g.uniform(0, 1, 260)
    curve = np.c_[4 * t - 2, 0.9 * np.sin(2.4 * np.pi * t)] + g.normal(0, 0.035, (260, 2))
    X = np.r_[curve, g.uniform([-2.3, -1.6], [2.3, 1.6], (420, 2))]
    return X, np.r_[np.ones(260, bool), np.zeros(420, bool)]


def test_preference_follows_the_paper(data):
    X, _ = data
    G = LAATGraph(X, 0.3)
    i = int(np.flatnonzero(G.valid)[10])
    a, b = G.indptr[i], G.indptr[i + 1]
    nb = X[G.col[a:b]]
    lam, V = np.linalg.eigh(np.cov(np.vstack([nb, X[i]]).T))
    d = nb - X[i]
    cos = np.abs(d @ V) / np.linalg.norm(d, axis=1, keepdims=True)
    E = (cos / cos.sum(1, keepdims=True)) @ (lam / lam.sum())
    assert np.allclose(G.E[a:b], E)
    assert np.allclose(G.Ebar[a:b], E / E.sum())
    assert lam.min() / lam.sum() - 1e-12 <= E.min() and E.max() <= lam.max() / lam.sum() + 1e-12


def test_probabilities_are_boltzmann_and_normalised(data):
    X, _ = data
    G = LAATGraph(X, 0.3)
    F = np.random.default_rng(1).uniform(0.5, 2, len(X))
    p = G.probabilities(F, 0.6, 40.0)
    assert np.allclose(np.add.reduceat(p, G.first[G.deg > 0]), 1)
    i = int(np.flatnonzero(G.valid)[3])
    a, b = G.indptr[i], G.indptr[i + 1]
    V = 0.4 * F[G.col[a:b]] / F[G.col[a:b]].sum() + 0.6 * G.Ebar[a:b]
    ref = np.exp(40.0 * (V - V.max()))
    assert np.allclose(p[a:b], ref / ref.sum())


@pytest.mark.parametrize("norm", ["sum", "max"])
def test_stationary_distribution(data, norm):
    X, _ = data
    G = LAATGraph(X, 0.3, norm=norm)
    pi = G.stationary(10.0)
    assert np.isclose(pi.sum(), 1) and (pi >= 0).all() and (pi[~G.valid] == 0).all()
    assert np.abs(G.P1.T @ pi - pi).max() < 1e-10
    assert G.n_components >= 1


def test_visits_follow_the_stationary_distribution(data):
    X, _ = data
    G = LAATGraph(X, 0.3)
    pi = G.stationary(30.0)
    tot = np.zeros(len(X))
    for s in range(60):
        tot += run(G, None, kappa=1.0, beta=30.0, n_ants=200, n_steps=20, n_epochs=1, zeta=0.0, phi=1.0, seed=s, start_p=pi) - 1.0
    assert np.corrcoef(tot, pi)[0, 1] > 0.98


def test_max_normalisation_does_not_depend_on_density():
    # two copies of the same neighbourhood shape, one with every neighbour duplicated nearby:
    # with the maximum the largest log-odds are the same, with the sum they shrink
    X = np.array([[0, 0], [1, 0], [-1, 0], [0, 0.3], [0, -0.3], [2.5, 0]], float)
    G = LAATGraph(X, 1.2, min_nb=2, norm="max")
    assert np.isclose(G.Ebar[G.indptr[0]:G.indptr[1]].max(), 1.0)
    Gs = LAATGraph(X, 1.2, min_nb=2, norm="sum")
    assert np.isclose(Gs.Ebar[Gs.indptr[0]:Gs.indptr[1]].sum(), 1.0)


def test_start_rules(data):
    X, _ = data
    G = LAATGraph(X, 0.3)
    for rule in laat_exact.START_RULES:
        w = G.start_distribution(rule, 10.0)
        assert np.isclose(w.sum(), 1) and (w >= 0).all() and (w[~G.valid] == 0).all()
    assert (G.start_distribution("median")[G.deg < np.median(G.deg)] == 0).all()
    with pytest.raises(ValueError):
        G.start_distribution("stationary")
    with pytest.raises(ValueError):
        G.start_distribution("nonsense")


def test_run_is_reproducible_and_finds_the_curve(data):
    X, truth = data
    G = LAATGraph(X, 0.3)
    f1 = run(G, "median", seed=3, n_epochs=5)
    f2 = run(G, "median", seed=3, n_epochs=5)
    assert np.array_equal(f1, f2)
    F = run(G, "median", kappa=0.8, beta=10.0, n_epochs=20, seed=0)
    assert F[truth].mean() > 2 * F[~truth].mean()


def test_walk_lengths():
    L = walk_lengths(np.array([1, 1, 1, 10, 100.0]), 500)
    assert L.min() >= 1 and abs(L.sum() - 500) <= 5 and L[-1] == L.max()
    G_len = lambda a: np.ones(len(a))
    assert np.array_equal(walk_lengths(G_len(np.arange(4)), 40), np.full(4, 10))


def test_individual_walk_lengths_keep_the_budget(data):
    X, _ = data
    G = LAATGraph(X, 0.3)
    pi = G.stationary(10.0)
    F = run(G, "median", beta=10.0, n_ants=50, n_steps=40, n_epochs=1, zeta=0.0, phi=1.0, seed=0,
            length=lambda a: 1.0 / np.maximum(pi[a], pi[pi > 0].min()))
    assert abs((F - 1).sum() - 50 * 40) <= 0.05 * 50 * 40
