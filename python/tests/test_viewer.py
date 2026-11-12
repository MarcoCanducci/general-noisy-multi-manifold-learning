import base64
import json
import re
import warnings

import numpy as np
import pytest

pytest.importorskip("plotly")

from examm import HierarchicalAGTM, ManifoldGraph, datasets  # noqa: E402
from examm.crawling import crawl  # noqa: E402
from examm.viewer import export_viewer, viewer_payload  # noqa: E402


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
