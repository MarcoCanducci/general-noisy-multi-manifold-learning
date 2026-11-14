import warnings

import numpy as np
import pytest

import manydream
from manydream import results


@pytest.fixture(scope="module")
def fitted():
    d = manydream.datasets.make_filaments(random_state=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        L = manydream.MultiManifoldLearner(r=0.1, dims=[1], random_state=0).fit(d.X)
    return d, L


def test_skeleton_lists_every_node_and_edge(fitted, tmp_path):
    d, L = fitted
    sk = results.skeleton(L)
    assert len(sk["nodes"]) == sum(m.graph.n_nodes for m in L.manifolds_)
    assert len(sk["edges"]) == sum(len(m.graph.edges) for m in L.manifolds_) and (sk["edges"]["length"] > 0).all()
    assert np.allclose(np.c_[sk["nodes"]["x1"], sk["nodes"]["x2"], sk["nodes"]["x3"]][:3], L.manifolds_[0].agtm.means_[:3])
    a, b = results.save_skeleton(L, tmp_path / "sk")
    assert a.read_text().startswith("structure,node,dim,x1,x2,x3") and b.read_text().startswith("structure,node_a,node_b,length")


def test_membership_separates_structure_from_background(fitted, tmp_path):
    d, L = fitted
    r = results.membership(L)
    P = r["proba"]
    assert P.shape == (len(d.X), len(L.manifolds_) + 1) and np.allclose(P.sum(1), 1, atol=1e-6)
    assert np.median(r["on_structure"][d.label > 0]) > 0.9 > 0.2 > np.median(r["on_structure"][d.label == 0])
    far = results.membership(L, d.X.max(0)[None] + 0.0)                 # a corner of the box: background
    assert far["structure"][0] == -1
    out = results.save_membership(L, tmp_path / "m.csv")
    assert out.read_text().splitlines()[0].startswith("row,structure,p_on_structure,p_background,p_structure_0")


def test_validate_counts_errors_and_matches_manifolds(fitted):
    d, L = fitted
    rep = results.validate(results.membership(L)["structure"], d.label, d.names)
    det = rep.detection
    assert det["tp"] + det["fn"] == (d.label > 0).sum() and det["fp"] + det["tn"] == (d.label == 0).sum()
    assert det["precision"] > 0.8 and det["recall"] > 0.8
    assert len(rep.manifolds) == 3 and all(m["recall"] > 0.6 for m in rep.manifolds)
    assert "false positives" in str(rep)
    # a hand-made case: one manifold split in two, one noise structure
    truth = np.r_[np.ones(100), np.zeros(50)].astype(int)
    pred = np.r_[np.zeros(60), np.ones(40), np.full(30, 2), np.full(20, -1)].astype(int)
    rep = results.validate(pred, truth)
    assert rep.detection == dict(tp=100, fp=30, fn=0, tn=20, precision=100 / 130, recall=1.0, f1=2 * (100 / 130) / (1 + 100 / 130))
    assert rep.manifolds[0]["pieces"] == 2 and rep.manifolds[0]["recall"] == 0.6 and rep.spurious[0]["structure"] == 2
