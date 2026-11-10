import warnings

import numpy as np
import pytest

import examm
from examm import MultiManifoldLearner, datasets


@pytest.fixture(scope="module")
def toroids():
    ds = datasets.make_three_toroids(scale=0.5, random_state=7)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        learner = MultiManifoldLearner(r=0.1, dims=[2], random_state=0, agtm_params={"max_iter": 15}).fit(ds.X)
    return ds, learner


def test_three_toroids_end_to_end(toroids):
    ds, L = toroids
    assert [m.dim for m in L.manifolds_] == [2, 2, 2]
    assert L.stale_stages == ()
    found = set()
    for k, m in enumerate(L.manifolds_):
        counts = np.bincount(ds.label[m.point_indices], minlength=4)
        assert counts[1:].max() / counts.sum() > 0.97
        found.add(int(counts[1:].argmax()))
        assert m.graph.degree.max() <= 4 and m.agtm is not None
    assert found == {0, 1, 2}
    attached = L.labels_ >= 0
    assert attached[ds.label > 0].mean() > 0.9 and attached[ds.label == 0].mean() < 0.2
    assert L.dim_index_.shape == (len(ds.X),) and set(np.unique(L.dim_index_)) <= {0, 1, 2, 3}
    # the density model assigns points to the right torus
    P = L.predict(ds.X[ds.label > 0][:3000])
    truth = ds.label[ds.label > 0][:3000]
    mapping = {k: int(np.bincount(ds.label[m.point_indices], minlength=4)[1:].argmax()) + 1
               for k, m in enumerate(L.manifolds_)}
    assert np.mean([mapping[p] == t for p, t in zip(P, truth)]) > 0.97
    assert np.isfinite(L.score(ds.X[:1000])) and "manifold 2" in L.summary()


def test_refit_recomputes_only_stale_stages(toroids):
    _, L = toroids
    pre, dim, graphs = L.preprocess_, L.dim_result_, L.graphs_
    L.set_params(epsilon=2)
    assert L.stale_stages == ("agtm",)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        L.run()
    assert L.preprocess_ is pre and L.dim_result_ is dim and L.graphs_ is graphs
    assert all(len(m.agtm.centres_) < m.graph.n_nodes / 3 for m in L.manifolds_)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        L.refit(beta=0.5)
    assert L.preprocess_ is pre and L.dim_result_ is dim and L.graphs_ is not graphs
    assert len(L.manifolds_) == 3
    L.set_params(r=0.11)
    assert L.stale_stages == examm.STAGES
    L.set_params(r=0.1)
    with pytest.raises(ValueError):
        L.set_params(nonsense=1)
    with pytest.raises(ValueError):
        L.refit(from_stage="nope")


def test_save_load_roundtrip(toroids, tmp_path):
    ds, L = toroids
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        L.run()
    examm.io.save(L, tmp_path / "model.examm")
    L2 = examm.io.load(tmp_path / "model.examm")
    assert np.allclose(L2.score_samples(ds.X[:500]), L.score_samples(ds.X[:500]))
    assert len(L2.manifolds_) == len(L.manifolds_)


def test_load_table_formats(tmp_path, rng):
    A = rng.normal(size=(20, 4))
    np.savetxt(tmp_path / "a.csv", A, delimiter=",", header="x,y,z,rho", comments="")
    np.savetxt(tmp_path / "b.txt", A)
    np.save(tmp_path / "c.npy", A)
    np.savez(tmp_path / "d.npz", data=A)
    B, names = examm.io.load_table(tmp_path / "a.csv")
    assert np.allclose(A, B) and names == ["x", "y", "z", "rho"]
    for f in ("b.txt", "c.npy", "d.npz"):
        B, names = examm.io.load_table(tmp_path / f)
        assert np.allclose(A, B) and names == ["x1", "x2", "x3", "x4"]


def test_other_configurations_and_errors(rng):
    X = datasets.circle(3000, 1.0, 2, rng) + 0.01 * rng.normal(size=(3000, 2))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        L = MultiManifoldLearner(r=0.15, preprocess=None, random_state=1).fit(X)
        assert len(L.manifolds_) == 1 and L.manifolds_[0].dim == 1
        E = MultiManifoldLearner(r=0.15, preprocess=None, graph="epsilon", attach="all", background=0.1,
                                 random_state=1).fit(X)
    assert len(E.manifolds_) >= 1 and E.model_.background == 0.1
    for bad in (dict(r=-1.0), dict(r=0.1, graph="x"), dict(r=0.1, attach="x"), dict(r=0.1, dims=[5])):
        with pytest.raises(ValueError):
            MultiManifoldLearner(**bad).fit(X)
    with pytest.raises(ValueError):
        MultiManifoldLearner(r=0.1).fit(np.full((10, 2), np.nan))
    with pytest.raises(RuntimeError):
        MultiManifoldLearner(r=0.1).run()
    # nothing to find: no crash, empty result
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        N = MultiManifoldLearner(r=0.01, preprocess=None).fit(rng.uniform(size=(500, 3)))
    assert N.manifolds_ == [] and N.model_ is None and (N.labels_ == -1).all()
    with pytest.raises(RuntimeError):
        N.predict(X)


def test_d1_recovers_manifolds_of_mixed_dimension():
    ds = datasets.make_d1(random_state=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        L = MultiManifoldLearner(r=0.1, random_state=0, agtm_params={"max_iter": 10}).fit(ds.X)
    owners = {}
    for m in L.manifolds_:
        counts = np.bincount(ds.label[m.point_indices], minlength=8)
        k = int(counts[1:].argmax()) + 1
        assert counts[k] / counts[1:].sum() > 0.95
        assert m.dim == {1: 1, 2: 1, 3: 2, 4: 2, 5: 2, 6: 2}.get(k, m.dim) or m.graph.n_nodes < 10
        owners.setdefault(k, []).append(m)
    # the four surfaces are each found as one 2-D manifold, the parabola as one curve
    for k in (3, 4, 5, 6):
        assert sum(m.dim == 2 for m in owners[k]) == 1
    assert sum(m.dim == 1 for m in owners[1]) == 1 and 1 <= len(owners[2]) <= 4
