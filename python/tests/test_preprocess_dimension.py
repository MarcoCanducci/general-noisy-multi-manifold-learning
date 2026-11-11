import numpy as np
import pytest

import examm
from examm import datasets
from examm.dimension import barycentric, dimensionality_index, fisher_distance
from examm.preprocess import SAF, DensityFilter, DensityThin, Pipeline, PreprocessResult, Preprocessor, resolve


def saf_reference(samples, original, h, mu, support):
    out = np.empty_like(samples)
    for i, v in enumerate(samples):
        d = np.linalg.norm(original - v, axis=1)
        m = d <= support
        w = np.exp(-d[m] ** 2 / (2 * h * h))
        avg = (w[:, None] * original[m]).sum(0) / w.sum() if w.sum() > 1e-20 else v
        ds = np.linalg.norm(samples - v, axis=1)
        m = (ds <= support) & (np.arange(len(samples)) != i)
        w = np.exp(-ds[m] ** 2 / (2 * h * h)) / np.maximum(ds[m], 1e-8)
        rep = (w[:, None] * (v - samples[m])).sum(0) / w.sum() if w.sum() > 1e-20 else 0
        out[i] = avg + mu * rep
    return out


@pytest.mark.parametrize("d", [2, 3, 5])
def test_saf_step_matches_reference(rng, d):
    X = rng.normal(size=(300, d))
    S = X + 0.01 * rng.normal(size=X.shape)
    got = examm._core.KDTree(X).saf_step(S, 0.4, 0.3, 0.8)
    assert np.allclose(got, saf_reference(S, X, 0.4, 0.3, 0.8))


def test_saf_reduces_transverse_noise(rng):
    t = rng.uniform(0, 2 * np.pi, 4000)
    X = np.column_stack([np.cos(t), np.sin(t), np.zeros_like(t)]) + 0.03 * rng.normal(size=(4000, 3))
    res = SAF().fit_transform(X, r=0.2)
    err = lambda P: np.abs(np.linalg.norm(P[:, :2], axis=1) - 1).mean() + np.abs(P[:, 2]).mean()  # noqa: E731
    assert err(res.diffused) < 0.5 * err(X)
    assert res.diffused.shape == X.shape and np.array_equal(res.indices, np.arange(len(X)))


def test_density_filter_is_order_independent_and_removes_noise(rng):
    ds = datasets.make_three_toroids(scale=0.5, noise_ratio=0.5, random_state=1)
    X = ds.X
    for mode in ("core", "core+border"):
        pipe = Pipeline([SAF(), DensityFilter(mode=mode)])
        a = pipe.fit_transform(X, r=0.1)
        perm = rng.permutation(len(X))
        b = pipe.fit_transform(X[perm], r=0.1)
        assert set(a.indices.tolist()) == set(perm[b.indices].tolist())
        assert np.allclose(a.noisy, X[a.indices])
        if mode == "core":
            core = a
    kept = np.zeros(len(X), bool)
    kept[core.indices] = True
    assert kept[ds.label > 0].mean() > 0.95     # manifold points survive
    assert kept[ds.label == 0].mean() < 0.25    # most background removed
    assert set(core.indices.tolist()) <= set(a.indices.tolist())


def test_custom_preprocessor_plugs_in(rng):
    class Halve(Preprocessor):
        def apply(self, state, r):
            return state.select(state.noisy[:, 0] > 0)

    X = rng.normal(size=(200, 3))
    out = Pipeline([Halve(), SAF(n_iter=1)]).fit_transform(X, r=1.0)
    assert (out.noisy[:, 0] > 0).all() and len(out) == (X[:, 0] > 0).sum()
    examm.preprocess.register("halve", Halve)
    assert isinstance(resolve("halve"), Halve) and "halve" in examm.preprocess.available()
    with pytest.raises(ValueError):
        resolve("does-not-exist")
    assert len(resolve(None).fit_transform(X, 1.0)) == len(X)


def test_barycentric_and_fisher():
    d = 5
    for j in range(1, d + 1):
        s = np.r_[np.full(j, 1 / j), np.zeros(d - j)]
        assert np.allclose(barycentric(s), np.eye(d)[j - 1])
    lam = np.array([0.5, 0.3, 0.2])
    a = barycentric(lam)
    assert np.isclose(a.sum(), 1) and np.allclose(a, [0.2, 0.2, 0.6])
    assert np.isclose(fisher_distance(np.eye(3)[0], np.eye(3)[1]), np.pi)
    assert fisher_distance(lam, lam) < 1e-7


@pytest.mark.parametrize("d", [3, 4, 6])
def test_dimension_index_on_flats(rng, d):
    for j in range(1, d + 1):
        X = np.zeros((6000, d))
        X[:, :j] = rng.uniform(-1, 1, (6000, j))
        X += 1e-3 * rng.normal(size=X.shape)
        # rotate so nothing is axis aligned
        Qm, _ = np.linalg.qr(rng.normal(size=(d, d)))
        r = 0.25 * (j ** 0.5) + 0.05 * j
        res = dimensionality_index(X @ Qm, r=r)
        inner = np.abs(X[:, :j]).max(1) < 0.6
        assert (res.index[inner] == j).mean() > 0.97, (d, j)
        assert (res.index_hard[inner] == j).mean() > 0.9, (d, j)


def test_dimension_index_marks_isolated_points_as_zero(rng):
    X = np.vstack([rng.uniform(size=(500, 3)) * 0.1, [[5.0, 5, 5]]])
    res = dimensionality_index(X, r=0.05)
    assert res.index[-1] == 0 and res.proba_smooth[-1].sum() == 0
    assert len(res.index) == len(X) and (res.index[:-1] > 0).mean() > 0.9


def test_d1_dimension_index_matches_paper_quality():
    ds = datasets.make_d1(random_state=0)
    sizes = np.bincount(ds.label)
    assert abs(sizes[0] - 2 * sizes[1:].sum()) <= 1 and 30000 < sizes[1:].sum() < 40000
    res = resolve("saf").fit_transform(ds.X, r=0.1)
    kept = np.isin(np.arange(len(ds.X)), res.indices)
    assert kept[ds.label == 0].mean() < 0.08          # paper: 5 % of the noise survives
    assert kept[ds.label > 0].mean() > 0.99
    di = dimensionality_index(res.diffused, 0.1)
    dim = ds.dim[res.indices]
    for j, want in ((1, 0.9), (2, 0.95), (3, 0.9)):   # paper Fig. 12: 0.98 / 0.99 / 0.95
        assert (di.index[dim == j] == j).mean() > want, j


def test_density_thinning_evens_out_the_sampling(rng):
    dense = rng.normal(size=(20000, 3)) * 0.2
    sparse = rng.uniform(-3, 3, size=(2000, 3))
    X = np.vstack([dense, sparse])
    out = DensityThin(target=20).fit_transform(X, r=0.3)
    tree = examm._core.KDTree(out.noisy)
    before = examm._core.KDTree(X).radius_count(X, 0.3)
    after = tree.radius_count(out.noisy, 0.3)
    assert before.max() > 1000 and np.quantile(after, 0.99) < 45        # dense core reduced to about the target
    kept = np.isin(np.arange(len(X)), out.indices)
    assert kept[before <= 20].all()                                      # sparse regions untouched
    assert np.allclose(out.noisy, X[out.indices]) and np.allclose(out.diffused, out.noisy)
    assert len(out.info["thin_probability"]) == len(out)
    again = DensityThin(target=20).fit_transform(X, r=0.3)
    assert np.array_equal(out.indices, again.indices)                    # reproducible by default
    other = DensityThin(target=20, random_state=1).fit_transform(X, r=0.3)
    assert not np.array_equal(out.indices, other.indices)
    piped = Pipeline([DensityThin(target=20), SAF(n_iter=2)]).fit_transform(X, r=0.3)
    assert np.array_equal(piped.indices, out.indices)
    assert isinstance(resolve("thin"), DensityThin)
    with pytest.raises(ValueError):
        DensityThin(target=0)
