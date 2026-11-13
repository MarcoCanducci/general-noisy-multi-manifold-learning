import base64
import json
import re
import warnings

import numpy as np
import pytest

pytest.importorskip("plotly")

from manydream import HierarchicalAGTM, ManifoldGraph, datasets  # noqa: E402
from manydream.crawling import crawl  # noqa: E402
from manydream.viewer import _graph_embedding, export_viewer, viewer_payload  # noqa: E402


def small_hierarchy(rng):
    x = rng.uniform(0, 4, 2000)
    lines = np.vstack([np.c_[x[:1000], np.full(1000, 0.2), np.zeros(1000)], np.c_[x[1000:], np.full(1000, -0.2), np.zeros(1000)]])
    c = datasets.circle(1000, 1.0, 3, rng) + [2, 3, 0]
    X = np.vstack([lines, c]) + 0.02 * rng.normal(size=(3000, 3))
    t = np.linspace(0, 4, 21)
    coarse = ManifoldGraph(np.c_[t, np.zeros(21), np.zeros(21)], [[i, i + 1] for i in range(20)], 1)
    fine = [ManifoldGraph(np.c_[t, np.full(21, s), np.zeros(21)], [[i, i + 1] for i in range(20)], 1) for s in (0.2, -0.2)]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        H = HierarchicalAGTM(background=0.01).set_data(X)
        p, _ = H.add_children([coarse, crawl(c, 0.2, 1, seed=0)])
        H.add_children(fine, p)
    return X, H


def test_payload_describes_tree_and_points(rng):
    X, H = small_hierarchy(rng)
    P = viewer_payload(X, {"a": H}, names=["u", "v", "w"])
    v = P["versions"][0]
    S = v["structures"]
    assert [s["kind"] for s in S] == ["group", "leaf", "leaf", "residual", "leaf"]
    assert [s["parent"] for s in S] == [-1, 0, 0, 0, -1] and S[3]["graph_of"] == 0
    sid = np.frombuffer(base64.b64decode(v["sid"]), np.int16)
    assert len(sid) == len(X) and sid.min() >= 1 and set(np.unique(sid)) <= {1, 2, 3, 4}
    assert [s["points"] for s in S] == [0] + [int((sid == k).sum()) for k in (1, 2, 3, 4)]
    assert (sid[:1000] == 1).mean() > 0.9 and (sid[1000:2000] == 2).mean() > 0.9 and (sid[2000:] == 4).mean() > 0.99
    for key in ("logdens", "dist", "d2", "pval", "rank", "bg"):
        assert len(np.frombuffer(base64.b64decode(v["scores"][key]), np.float32)) == len(X)
    nodes = np.frombuffer(base64.b64decode(S[4]["nodes"]), np.float32).reshape(-1, 3)
    assert len(nodes) == S[4]["n_nodes"] and np.allclose(np.frombuffer(base64.b64decode(P["X"]), np.float32).reshape(-1, 3), X, atol=1e-5)
    with pytest.raises(ValueError):
        viewer_payload(X, {})
    with pytest.raises(ValueError):
        viewer_payload(X, {"a": H}, names=["only-one"])


def test_export_writes_a_self_contained_page(rng, tmp_path):
    X, H = small_hierarchy(rng)
    out = export_viewer(tmp_path / "v.html", X, H, title="my viewer", plotlyjs="cdn")
    html = out.read_text()
    assert "<title>my viewer</title>" in html and "cdn.plot.ly" in html and "/*__DATA__*/" not in html
    data = json.loads(re.search(r"const DATA = (\{.*?\});\nconst PAL", html, re.S).group(1))
    assert data["n"] == len(X) and data["versions"][0]["name"] == "model"
    inline = export_viewer(tmp_path / "w.html", X, {"a": H, "b": H})
    assert inline.stat().st_size > 3_000_000 and "cdn.plot.ly/plotly-" not in inline.read_text()[:2000]
    with pytest.raises(ValueError):
        export_viewer(tmp_path / "x.html", X, H, plotlyjs="nope")


def test_extra_labels_are_embedded(rng, tmp_path):
    X, H = small_hierarchy(rng)
    z = rng.uniform(size=len(X))
    cls = np.where(X[:, 1] > 1, "ring", "lines")
    P = viewer_payload(X, {"a": H}, extra={"redshift": z, "class": cls})
    num, cat = P["extra"]
    assert num["kind"] == "num" and np.allclose(np.frombuffer(base64.b64decode(num["values"]), np.float32), z, atol=1e-6)
    assert cat["kind"] == "cat" and cat["cats"] == ["lines", "ring"]
    codes = np.frombuffer(base64.b64decode(cat["values"]), np.int32)
    assert np.array_equal(np.array(cat["cats"])[codes], cls)
    out = export_viewer(tmp_path / "e.html", X, H, plotlyjs="cdn", extra={"redshift": z})
    assert '"name":"redshift"' in out.read_text() and 'id="xfile"' in out.read_text()
    with pytest.raises(ValueError):
        viewer_payload(X, {"a": H}, extra={"bad": z[:10]})


def test_graph_embedding_unfolds_structures_in_their_own_dimension(rng):
    # a curved surface: embedding distances follow the distances measured along the surface
    u, v = rng.uniform(0, 2, (2, 6000))
    g = crawl(np.c_[u, v, np.sin(2 * u)], 0.15, 2, seed=0)
    emb, arc = _graph_embedding(g.nodes, g.edges, 2)
    assert emb.shape == (g.n_nodes, 2) and len(arc) == 0
    t = np.linspace(0, 2, 400)
    s_of_u = np.interp(g.nodes[:, 0], t, np.r_[0, np.cumsum(np.hypot(np.diff(t), np.diff(np.sin(2 * t))))])
    true = np.c_[s_of_u, g.nodes[:, 1]]
    i, j = rng.integers(0, g.n_nodes, (2, 2000))
    d_emb, d_true = np.linalg.norm(emb[i] - emb[j], axis=1), np.linalg.norm(true[i] - true[j], axis=1)
    assert np.corrcoef(d_emb, d_true)[0, 1] > 0.95
    lm, _ = _graph_embedding(g.nodes, g.edges, 2, max_landmarks=200)       # landmark variant agrees
    assert np.corrcoef(np.linalg.norm(lm[i] - lm[j], axis=1), d_true)[0, 1] > 0.95
    # an open curve gets an arc length; 3-D and 4-D structures get 3 and 4 coordinates
    line = ManifoldGraph(np.c_[np.linspace(0, 3, 16), np.zeros(16), np.zeros(16)], [[k, k + 1] for k in range(15)], 1)
    e1, a1 = _graph_embedding(line.nodes, line.edges, 1)
    assert np.allclose(np.sort(a1), np.linspace(0, 3, 16)) and e1.shape == (16, 2)
    assert _graph_embedding(g.nodes, g.edges, 3)[0].shape[1] == 3 and _graph_embedding(g.nodes, g.edges, 6)[0].shape[1] == 4
    X, H = small_hierarchy(rng)
    S = viewer_payload(X, {"a": H})["versions"][0]["structures"]
    assert S[1]["m"] == 2 and all(k in S[1] for k in ("emb", "arc", "curv", "zeta", "degen", "noise"))
    assert len(np.frombuffer(base64.b64decode(S[1]["curv"]), np.float32)) == 20
    # covariances for the density surfaces: upper triangles, row by row
    tri = np.frombuffer(base64.b64decode(S[1]["cov"]), np.float32).reshape(S[1]["n_nodes"], 6)
    leaf = H.root_.children[0].children[0].model
    assert np.allclose(tri, leaf.covariances_[:, np.triu_indices(3)[0], np.triu_indices(3)[1]], rtol=1e-5, atol=1e-8)
    assert S[3]["cov"] == "" and S[3]["graph_of"] == 0


def test_points_fill_the_layout_of_a_surface(rng):
    from manydream.viewer import _point_layout

    # a flat lattice in the plane of a 5-D space, laid out by the identity: points must land on themselves
    gx, gy = np.meshgrid(np.arange(6.0), np.arange(6.0), indexing="ij")
    emb = np.c_[gx.ravel(), gy.ravel()]
    means = np.c_[emb, np.zeros((36, 3))]
    idx = np.arange(36).reshape(6, 6)
    edges = np.r_[np.c_[idx[:-1].ravel(), idx[1:].ravel()], np.c_[idx[:, :-1].ravel(), idx[:, 1:].ravel()]]
    P = np.c_[rng.uniform(0, 5, (200, 2)), rng.normal(0, 0.3, (200, 3))]     # off-plane noise is dropped
    node = np.argmin(((P[:, None, :2] - emb[None]) ** 2).sum(-1), axis=1)
    out = _point_layout(P, means, emb, edges, node, 2)
    assert np.allclose(out, P[:, :2], atol=1e-8)


def test_cloud_embedding_recovers_a_flat_sheet(rng):
    from manydream.viewer import _cloud_embedding

    P = np.c_[rng.uniform(0, 4, (1500, 2)), np.zeros((1500, 3))]
    means = np.c_[np.array([[1.0, 1.0], [3.0, 1.0], [1.0, 3.0], [3.0, 3.0]]), np.zeros((4, 3))]
    e, p = _cloud_embedding(means, np.array([[0, 1], [2, 3]]), P, 2, max_landmarks=300)
    assert e.shape == (4, 2) and p.shape == (1500, 2)
    d_true = np.linalg.norm(P[:200, None, :2] - P[None, :200, :2], axis=-1)
    d_emb = np.linalg.norm(p[:200, None] - p[None, :200], axis=-1)
    assert np.corrcoef(d_true.ravel(), d_emb.ravel())[0, 1] > 0.98


def test_dimension_labels_follow_the_detection_runs(rng):
    import manydream
    from manydream.viewer import _dimension_labels, viewer_payload

    X = np.r_[np.c_[rng.uniform(0, 3, (900, 2)), rng.normal(0, 0.01, 900)], rng.uniform(-1, 4, (60, 3))]
    H = manydream.build_hierarchy(X, [0.3], dims=[[1, 2]], random_state=0)
    hard, smooth = _dimension_labels(H, len(X))
    L = H.detections_[0]["learner"]
    assert np.array_equal(smooth[L.preprocess_.indices], L.dim_result_.index)
    assert np.array_equal(hard[L.preprocess_.indices], L.dim_result_.index_hard * (L.dim_result_.index > 0))
    removed = np.setdiff1d(np.arange(len(X)), L.preprocess_.indices)
    assert not hard[removed].any() and not smooth[removed].any()
    v = viewer_payload(X, {"a": H})["versions"][0]
    assert v["idxh"] and v["idxs"]
    del H.detections_
    assert _dimension_labels(H, len(X)) is None


def test_point_status_accounts_for_every_point(rng):
    import manydream
    from manydream.viewer import _point_status, viewer_payload

    X = np.r_[np.c_[rng.uniform(0, 3, (2500, 2)), rng.normal(0, 0.01, 2500)], rng.uniform(-1, 4, (80, 3))]
    pre = manydream.Pipeline([manydream.DensityThin(radius=0.2, target=20, random_state=0), manydream.SAF(), manydream.DensityFilter()])
    H = manydream.build_hierarchy(X, [0.3], preprocess=[pre], dims=[[1, 2]], random_state=0)
    L = H.detections_[0]["learner"]
    assert [n for n, _ in L.preprocess_.info["steps"]] == ["DensityThin", "SAF", "DensityFilter"]
    status, names, text = _point_status(H, len(X))
    assert status.min() >= 0 and len(status) == len(X)
    removed = np.setdiff1d(np.arange(len(X)), L.preprocess_.indices)
    assert (status[removed] > 0).all()
    thin = names.index("thinned out (density subsampling)")
    assert (status == thin).sum() == len(X) - len(L.preprocess_.info["steps"][0][1])
    assert "used: on a structure's graph" in text
    v = viewer_payload(X, {"a": H})["versions"][0]
    assert v["status"] and v["status_names"] == names
