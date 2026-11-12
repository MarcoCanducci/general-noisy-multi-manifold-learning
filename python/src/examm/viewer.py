"""Standalone interactive viewer (one HTML file) for one or several fitted hierarchies.

    examm.viewer.export_viewer("viewer.html", X, {"run A": H_a, "run B": H_b}, names=feature_names)

The page needs no server: open it in a browser. It shows the points in 3-D (rotate, zoom, pick
the three features), the graphs of the structures, a tree to switch structures, parents or leaves
on and off, and an outlier selector working on any of the per-point scores, by percentile or by
an explicit threshold, over all shown points or within each structure.
"""
from __future__ import annotations

import base64
import json
from pathlib import Path

import numpy as np

from ._util import as_points

__all__ = ["export_viewer", "viewer_payload"]

SCORES = [
    # key, label, low values are outliers?
    ("logdens", "log density (global model)", True),
    ("dist", "orthogonal distance from the structure", False),
    ("d2", "orthogonal Mahalanobis distance (squared)", False),
    ("pval", "model tail probability of the orthogonal offset", True),
    ("rank", "regional rank of the orthogonal distance", False),
    ("bg", "background probability", False),
]


def _b64(a, dtype) -> str:
    return base64.b64encode(np.ascontiguousarray(a, dtype=dtype).tobytes()).decode("ascii")


def _graph_embedding(means: np.ndarray, edges: np.ndarray, dim: int, max_landmarks: int = 1500, seed: int = 0):
    """Intrinsic coordinates of a graph, for drawing it in its own dimension.

    The graph is unfolded by classical multidimensional scaling of the geodesic distances
    measured along its edges (Isomap): ``m = min(max(dim, 2), 4)`` coordinates, so a surface is
    drawn in the plane, a 3-D structure in space, and a 4-D one through 3-D projections of four
    coordinates. Large graphs use landmark MDS. For chains the arc length of every node is
    returned as well (empty array otherwise).
    """
    from scipy.linalg import eigh
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import shortest_path

    n = len(means)
    m = int(min(max(dim, 2), 4))
    edges = np.asarray(edges, dtype=np.int64).reshape(-1, 2)
    arc = np.zeros(0)
    if dim == 1 and n >= 2:
        try:
            from .crawling import ManifoldGraph
            from .viz import chain_order

            order, s_ = chain_order(ManifoldGraph(means, edges, 1))
            if len(order) == n:
                arc = np.zeros(n)
                arc[order] = s_
        except Exception:  # noqa: BLE001 - the arc is optional
            arc = np.zeros(0)
    if n <= m or len(edges) == 0:
        return np.zeros((n, m)), arc
    w = np.linalg.norm(means[edges[:, 0]] - means[edges[:, 1]], axis=1) + 1e-12
    A = coo_matrix((np.r_[w, w], (np.r_[edges[:, 0], edges[:, 1]], np.r_[edges[:, 1], edges[:, 0]])), shape=(n, n)).tocsr()
    land = np.arange(n) if n <= max_landmarks else np.sort(np.random.default_rng(seed).choice(n, max_landmarks, replace=False))
    Dl = shortest_path(A, method="D", directed=False, indices=land)          # (L, n)
    finite = np.isfinite(Dl)
    Dl[~finite] = 1.5 * Dl[finite].max() if finite.any() else 1.0            # disconnected pieces: push apart
    D2 = Dl[:, land] ** 2
    L = len(land)
    J = np.eye(L) - 1.0 / L
    vals, vecs = eigh(-0.5 * J @ D2 @ J, subset_by_index=[max(L - m, 0), L - 1])
    vals, vecs = vals[::-1], vecs[:, ::-1]
    pos = np.maximum(vals, 1e-12)
    if n == L:
        emb = vecs * np.sqrt(pos)
    else:                                                                   # landmark MDS triangulation
        emb = (-0.5 * (vecs / np.sqrt(pos)).T @ (Dl ** 2 - D2.mean(1, keepdims=True))).T
    out = np.zeros((n, m))
    out[:, :emb.shape[1]] = emb
    return out, arc


def _point_layout(X, means, emb, edges, node, dim: int) -> np.ndarray:
    """Positions of data points in the intrinsic layout of their graph.

    Around every node the layout is locally linear: the offsets to the nodes within two hops,
    restricted to their ``dim`` leading directions (the local tangent space), are mapped on the
    corresponding offsets in the layout. A point is drawn at its node plus the image of its own
    offset, so points fill the faces of the graph and the component orthogonal to the structure
    is dropped.
    """
    n, m = emb.shape
    out = emb[node].astype(float)
    edges = np.asarray(edges, dtype=np.int64).reshape(-1, 2)
    if n < 2 or len(edges) == 0 or len(node) == 0:
        return out
    nbrs = [[] for _ in range(n)]
    for a, b in edges:
        nbrs[a].append(b)
        nbrs[b].append(a)
    order = np.argsort(node, kind="stable")
    starts = np.flatnonzero(np.r_[True, np.diff(node[order]) != 0, True])
    for s0, s1 in zip(starts[:-1], starts[1:]):
        rows = order[s0:s1]
        a = int(node[rows[0]])
        ring = set(nbrs[a])
        for b in nbrs[a]:
            ring.update(nbrs[b])
        ring.discard(a)
        if not ring:
            continue
        nb = np.fromiter(ring, dtype=np.int64)
        A, B = means[nb] - means[a], emb[nb] - emb[a]
        U, S, Vt = np.linalg.svd(A, full_matrices=False)
        k = int(min(dim, np.sum(S > 1e-8 * S[0]))) if S[0] > 0 else 0
        if k == 0:
            continue
        L = (Vt[:k].T / S[:k]) @ (U[:, :k].T @ B)                 # (D, m)
        off = (X[rows] - means[a]) @ L
        lim = np.linalg.norm(B, axis=1).max()                     # stay within the node's neighbourhood
        norm = np.linalg.norm(off, axis=1)
        off *= np.minimum(1.0, lim / np.maximum(norm, 1e-300))[:, None]
        out[rows] += off
    return out


def _version_payload(H, X, hops: int, min_region: int) -> dict:
    res = H.regional_deviation(X, hops=hops, min_region=min_region)
    leaves = H.leaves()
    structures, node_id = [], {}
    iu = np.triu_indices(X.shape[1])            # covariances are shipped as upper triangles, row by row

    def walk(node, parent_id):
        for c in node.children:
            sid = len(structures)
            node_id[id(c)] = sid
            g = c.model.graph
            emb, arc = _graph_embedding(c.model.means_, g.edges, int(g.dim))
            cov = c.model.covariances_
            degen = np.zeros(g.n_nodes, dtype=np.float32)
            degen[c.model.diagnostics_["degenerate"]] = 1.0
            structures.append({"id": sid, "label": c.label, "parent": parent_id, "level": c.level, "dim": int(g.dim),
                               "kind": "group" if c.children else "leaf", "weight": float(c.weight),
                               "nodes": _b64(c.model.means_, np.float32), "n_nodes": int(g.n_nodes),
                               "edges": _b64(g.edges, np.int32), "graph_of": sid,
                               "emb": _b64(emb, np.float32), "m": int(emb.shape[1]),
                               "arc": _b64(arc, np.float32) if len(arc) else "",
                               "curv": _b64(c.model.edge_curvature(), np.float32),
                               "zeta": _b64(c.model.zeta_, np.float32), "degen": _b64(degen, np.float32),
                               "noise": _b64(np.sqrt(np.einsum("kaa->k", cov) / cov.shape[1]), np.float32),
                               "cov": _b64(cov[:, iu[0], iu[1]], np.float32)})
            walk(c, sid)
            if c.children and c.pi_residual > 0:
                rid = len(structures)
                structures.append({"id": rid, "label": f"{c.label} (residual)", "parent": sid, "level": c.level + 1,
                                   "dim": int(g.dim), "kind": "residual", "weight": float(c.weight * c.pi_residual),
                                   "nodes": "", "n_nodes": 0, "edges": "", "graph_of": sid,
                                   "emb": "", "m": 0, "arc": "", "curv": "", "zeta": "", "degen": "", "noise": "", "cov": ""})
                node_id[("res", id(c))] = rid

    walk(H.root_, -1)
    # map leaf index (into H.leaves()) to structure id; background leaves never win in regional_deviation
    lookup = np.full(len(leaves), -1)
    for k, (n, kind, _) in enumerate(leaves):
        if kind != "background":
            lookup[k] = node_id[("res", id(n))] if kind == "residual" else node_id[id(n)]
    sid = lookup[res["leaf"]]
    # where every point is drawn on the graph page: up to four intrinsic coordinates
    ppos = np.zeros((len(X), 4), dtype=np.float32)
    graph_of = np.array([s["graph_of"] for s in structures] + [-1])
    gid = graph_of[sid]
    models = {node_id[k]: None for k in node_id if not isinstance(k, tuple)}

    def collect(node):
        for c in node.children:
            models[node_id[id(c)]] = c.model
            collect(c)

    collect(H.root_)
    for g_, model in models.items():
        rows = np.flatnonzero(gid == g_)
        if len(rows) == 0:
            continue
        emb = _graph_embedding(model.means_, model.graph.edges, int(model.graph.dim))[0]
        ppos[rows, :emb.shape[1]] = _point_layout(X[rows], model.means_, emb, model.graph.edges,
                                                  np.asarray(res["node"])[rows], int(model.graph.dim))
    counts = np.bincount(sid[sid >= 0], minlength=len(structures))
    for s in structures:
        s["points"] = int(counts[s["id"]])
    return {
        "structures": structures,
        "sid": _b64(sid, np.int16),
        "node": _b64(res["node"], np.int32),
        "toward": _b64(res["toward"], np.int32),
        "frac": _b64(res["fraction"], np.float32),
        "ppos": _b64(ppos, np.float32),
        "scores": {"logdens": _b64(H.score_samples(X), np.float32), "dist": _b64(res["distance"], np.float32),
                   "d2": _b64(res["d2"], np.float32), "pval": _b64(res["p_value"], np.float32),
                   "rank": _b64(res["regional_rank"], np.float32), "bg": _b64(res["background_probability"], np.float32)},
        "summary": H.summary(),
    }


def _extra_payload(extra, n: int) -> list:
    out = []
    for name, values in (extra or {}).items():
        a = np.asarray(values)
        if a.shape != (n,):
            raise ValueError(f"extra column {name!r} must hold one value per point")
        if np.issubdtype(a.dtype, np.number) or a.dtype == bool:
            out.append({"name": str(name), "kind": "num", "values": _b64(a.astype(float), np.float32), "cats": []})
        else:
            cats, codes = np.unique(a.astype(str), return_inverse=True)
            out.append({"name": str(name), "kind": "cat", "values": _b64(codes, np.int32), "cats": cats.tolist()})
    return out


def viewer_payload(X, versions: dict, names=None, hops: int = 2, min_region: int = 30, extra=None) -> dict:
    """The data embedded in the viewer page (JSON-serialisable)."""
    X = as_points(X)
    n, d = X.shape
    names = list(names) if names is not None else [f"x{k + 1}" for k in range(d)]
    if len(names) != d:
        raise ValueError("names must have one entry per feature")
    if not versions:
        raise ValueError("at least one version is required")
    return {"n": n, "d": d, "names": names, "X": _b64(X, np.float32), "extra": _extra_payload(extra, n),
            "scores": [{"key": k, "label": l, "low": low} for k, l, low in SCORES],
            "versions": [{"name": str(name), **_version_payload(H, X, hops, min_region)} for name, H in versions.items()]}


def export_viewer(path, X, versions, names=None, title: str = "examm viewer", plotlyjs: str = "inline",
                  hops: int = 2, min_region: int = 30, extra=None) -> Path:
    """Write the standalone viewer page.

    Parameters
    ----------
    path : file to write (``.html``)
    X : (n, d) data
    versions : dict name -> fitted :class:`examm.HierarchicalAGTM` (or a single hierarchy)
    names : feature names
    plotlyjs : ``"inline"`` embeds plotly.js in the page (works offline, about 5 MB);
        ``"cdn"`` loads it from the plotly CDN.
    hops, min_region : passed to :meth:`HierarchicalAGTM.regional_deviation`
    extra : dict name -> (n,) array, optional
        Additional per-point labels (numeric, or text/categorical) offered for colouring the
        points and the graph nodes, e.g. physical parameters. More can be added later in the
        page itself by loading a CSV file.
    """
    if not isinstance(versions, dict):
        versions = {"model": versions}
    payload = viewer_payload(X, versions, names, hops, min_region, extra)
    if plotlyjs == "inline":
        from plotly.offline import get_plotlyjs

        lib = "<script>" + get_plotlyjs() + "</script>"
    elif plotlyjs == "cdn":
        lib = '<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>'
    else:
        raise ValueError("plotlyjs must be 'inline' or 'cdn'")
    html = (_TEMPLATE.replace("/*__DATA__*/", json.dumps(payload, separators=(",", ":")))
            .replace("<!--__PLOTLY__-->", lib).replace("__TITLE__", title))
    path = Path(path)
    path.write_text(html, encoding="utf-8")
    return path


_TEMPLATE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
:root{--ink:#0b0b0b;--muted:#52514e;--line:#dddcd6;--bg:#fcfcfb;--panel:#f4f3ef;--accent:#2a78d6;--out:#e34948}
*{box-sizing:border-box}html,body{height:100%;margin:0}
body{font:13px/1.4 system-ui,-apple-system,Segoe UI,Roboto,sans-serif;color:var(--ink);background:var(--bg);display:flex}
#side{width:350px;min-width:350px;height:100%;overflow:auto;background:var(--panel);border-right:1px solid var(--line);padding:12px 14px}
#main{flex:1;min-width:0;height:100%;display:flex;flex-direction:column;overflow:hidden;position:relative}
#loc{position:absolute;left:12px;bottom:42px;width:340px;height:300px;z-index:3;background:#ffffff;border:1px solid var(--line);border-radius:6px;box-shadow:0 2px 10px rgba(0,0,0,.12);overflow:hidden}
#loc.big{width:620px;height:520px}#loc.off{display:none}
#lochead{height:24px;line-height:24px;padding:0 8px;font-size:11px;color:var(--muted);border-bottom:1px solid var(--line);display:flex;justify-content:space-between;cursor:pointer;user-select:none}
#locplot{position:absolute;top:25px;left:0;right:0;bottom:0}
#plot{flex:1 1 0;min-height:0;height:0;overflow:hidden;position:relative}
#status{flex:none;height:30px;line-height:18px;background:var(--bg);position:relative;z-index:2;padding:6px 12px;border-top:1px solid var(--line);color:var(--muted);font-size:12px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
h1{font-size:15px;margin:0 0 8px}h2{font-size:11px;text-transform:uppercase;letter-spacing:.06em;color:var(--muted);margin:16px 0 6px}
label{display:block;margin:6px 0 2px;color:var(--muted);font-size:12px}
select,input[type=number]{width:100%;padding:4px 6px;border:1px solid var(--line);border-radius:4px;background:#fff;font:inherit}
.row{display:flex;gap:6px;align-items:center}.row>*{flex:1}.row.tight>*{flex:none}
button{font:inherit;padding:4px 8px;border:1px solid var(--line);border-radius:4px;background:#fff;cursor:pointer}
button:hover{border-color:var(--accent)}
#tabs{flex:none;display:flex;gap:4px;padding:8px 12px 0;border-bottom:1px solid var(--line);background:var(--panel)}
.tab{border:1px solid var(--line);border-bottom:none;border-radius:6px 6px 0 0;background:#e9e8e3;padding:6px 14px}
.tab.active{background:var(--bg);font-weight:600;position:relative;top:1px}
#dimtog{display:flex;flex-wrap:wrap;gap:10px;margin:0 0 6px}
#dimtog label{display:flex;align-items:center;gap:4px;margin:0;color:var(--ink)}
#catlegend{display:flex;flex-wrap:wrap;gap:4px 10px;margin-top:4px;font-size:11px}#catlegend span{display:flex;align-items:center;gap:4px}
.viewonly{display:none}body.v3d .only3d{display:block}body.vgraph .onlygraph{display:block}
#tree{border:1px solid var(--line);border-radius:4px;background:#fff;max-height:300px;overflow:auto;padding:4px 0}
.node{display:flex;align-items:center;gap:5px;padding:1px 6px;white-space:nowrap;cursor:default}
.node:hover{background:#f1f5fb}.node input{margin:0}
.sw{width:10px;height:10px;border-radius:2px;flex:none}.cnt{color:var(--muted);font-size:11px;margin-left:auto;padding-left:8px}
.kind{color:var(--muted);font-size:11px}
.check{display:flex;align-items:center;gap:6px;margin:4px 0;color:var(--ink)}.check input{margin:0}
#outinfo{margin-top:6px;padding:6px 8px;background:#fff;border:1px solid var(--line);border-radius:4px;font-size:12px}
#outinfo b{color:var(--out)}
details{margin-top:10px}summary{cursor:pointer;color:var(--muted);font-size:12px}
pre{font-size:10.5px;white-space:pre;overflow:auto;max-height:240px;background:#fff;border:1px solid var(--line);border-radius:4px;padding:6px}
.hint{color:var(--muted);font-size:11px;margin:4px 0 0}
</style><!--__PLOTLY__--></head>
<body>
<div id="side">
  <h1>__TITLE__</h1>
  <label for="version">Version</label><select id="version"></select>
  <h2>View</h2>
  <div class="row"><div><label for="ax0">x</label><select id="ax0"></select></div><div><label for="ax1">y</label><select id="ax1"></select></div><div><label for="ax2">z</label><select id="ax2"></select></div></div>
  <label for="colour">Colour points by</label><select id="colour"></select>
  <div class="row"><div><label for="psize">Point size</label><input id="psize" type="range" min="1" max="6" step="0.5" value="2"></div>
  <div><label for="popac">Opacity</label><input id="popac" type="range" min="0.1" max="1" step="0.05" value="0.7"></div></div>
  <label class="check"><input type="checkbox" id="showpts" checked> points</label>
  <label class="check"><input type="checkbox" id="showgraphs" checked> graphs of the selected structures</label>
  <div class="viewonly onlygraph">
  <h2>Graph view</h2>
  <label for="gstruct">Structure</label><select id="gstruct"></select>
  <label for="gcolour">Colour nodes by</label><select id="gcolour"></select>
  <label for="gedge">Colour edges by</label><select id="gedge"><option value="none">nothing (grey)</option><option value="curv">curvature along the edge</option><option value="len">edge length in data space</option><option value="stretch">distortion of the layout</option></select>
  <div id="gprojwrap" style="display:none"><label for="gproj">Projection (three of the four intrinsic coordinates)</label><select id="gproj"><option value="123">coordinates 1, 2, 3</option><option value="124">coordinates 1, 2, 4</option><option value="134">coordinates 1, 3, 4</option><option value="234">coordinates 2, 3, 4</option></select></div>
  <label class="check"><input type="checkbox" id="gsize" checked> node size by number of points</label>
  <label for="gy">Vertical axis (curves)</label><select id="gy"></select>
  <label class="check"><input type="checkbox" id="gpts" checked> points of the structure</label>
  <label class="check"><input type="checkbox" id="gfill" checked> fill the structure (points at their own position)</label>
  <label class="check"><input type="checkbox" id="gloc" checked> locator: where the structure sits in the data</label>
  <p class="hint">Each structure is drawn in its own dimension, unfolded along its graph (geodesic distances along the edges): a curve as a line with a profile above it, a surface in the plane, a 3-D structure in space, a 4-D one through 3-D projections. Points are drawn at their own position within the structure (their offset from the nearest node, along the structure, carried into the layout); untick "fill the structure" to snap them onto the nearest edge instead.</p>
  </div>
  <h2>Structures</h2>
  <div id="dimtog"></div>
  <div class="row tight" style="margin-bottom:6px"><button id="all">All</button><button id="none">None</button><button id="leavesonly">Leaves only</button><button id="top">Residuals only</button></div>
  <div id="tree"></div>
  <p class="hint">Ticking a parent ticks everything under it. "(residual)" holds the points of a parent that none of its children explains; its graph is the parent's.</p>
  <div class="viewonly only3d">
  <h2>Model density</h2>
  <label for="isomode">Iso-surfaces (iso-curves in 2-D) of the fitted model</label>
  <select id="isomode"><option value="off">off</option><option value="all">the ticked structures together (all ticked = whole hierarchy)</option><option value="each">each ticked structure separately</option></select>
  <div class="row tight" style="flex-wrap:wrap;gap:4px 12px;margin-top:4px">
    <label class="check"><input type="checkbox" class="isolvl" value="0.5"> 50 %</label><label class="check"><input type="checkbox" class="isolvl" value="0.9" checked> 90 %</label>
    <label class="check"><input type="checkbox" class="isolvl" value="0.99" checked> 99 %</label><label class="check"><input type="checkbox" class="isolvl" value="0.999"> 99.9 %</label></div>
  <label for="isores">Resolution</label><select id="isores"><option value="32">coarse</option><option value="48" selected>medium</option><option value="64">fine</option></select>
  <p class="hint">Each surface encloses the stated share of the model's probability mass, for the model projected onto the shown features (its exact marginal). Points outside the outer surface are unlikely under the model <i>in this projection</i>; the outlier scores use all features, so a point can be an outlier without leaving the surfaces here. Choose "(none)" for z to get a 2-D view with iso-curves.</p>
  </div>
  <h2>Extra labels</h2>
  <input type="file" id="xfile" accept=".csv,.tsv,.txt">
  <p class="hint">Load a CSV with a header line to colour by your own per-point values (numbers or categories): either one line per point in the original order, or with a column named <b>row</b> holding the point index. The file stays in your browser.</p>
  <div id="xinfo" class="hint"></div><div id="catlegend"></div>
  <h2>Outliers</h2>
  <label for="score">Score</label><select id="score"></select>
  <label class="check"><input type="radio" name="mode" value="pct" checked> by percentile</label>
  <div class="row"><input id="pct" type="range" min="0.05" max="20" step="0.05" value="1"><input id="pctnum" type="number" min="0.01" max="50" step="0.05" value="1" style="max-width:90px"></div>
  <label class="check"><input type="radio" name="mode" value="thr"> by threshold</label>
  <div class="row"><input id="thr" type="range" step="any"><input id="thrnum" type="number" step="any" style="max-width:110px"></div>
  <label for="scope">Applied</label><select id="scope"><option value="shown">over all shown points</option><option value="each">within each shown structure</option></select>
  <label class="check"><input type="checkbox" id="hl" checked> highlight outliers</label>
  <label class="check"><input type="checkbox" id="only"> show outliers only</label>
  <div id="outinfo"></div>
  <div class="row tight" style="margin-top:6px"><button id="csv">Download outliers (CSV)</button><button id="csvall">Download shown points (CSV)</button></div>
  <details><summary>Model summary of this version</summary><pre id="summary"></pre></details>
</div>
<div id="main"><div id="tabs"><button class="tab active" data-view="3d">3-D view</button><button class="tab" data-view="graph">Graph view (one structure, unfolded)</button></div>
<div id="plot"></div><div id="status"></div>
<div id="loc" class="viewonly onlygraph"><div id="lochead" title="click to enlarge or shrink"><span id="loctitle"></span><span>&#x2922;</span></div><div id="locplot"></div></div></div>
<script>
const DATA = /*__DATA__*/;
const PAL = ["#2a78d6","#eb6834","#1baf7a","#eda100","#e87ba4","#008300","#6250d6","#e34948"];
const DIMCOL = ["#b9b8b2","#e34948","#2a78d6","#1baf7a","#eda100","#6250d6","#e87ba4","#008300","#52514e"];
function dec(b64, T){ if(!b64) return new T(0); const s=atob(b64), u=new Uint8Array(s.length); for(let i=0;i<s.length;i++)u[i]=s.charCodeAt(i); return new T(u.buffer); }
function colourOf(i){ if(i<PAL.length) return PAL[i]; const h=(i*137.508)%360; return `hsl(${h.toFixed(0)},62%,${i%2?44:56}%)`; }
const N=DATA.n, D=DATA.d, X=dec(DATA.X, Float32Array);
const V = DATA.versions.map(v => ({name:v.name, summary:v.summary, sid:dec(v.sid,Int16Array), node:dec(v.node,Int32Array),
  toward:dec(v.toward,Int32Array), frac:dec(v.frac,Float32Array), ppos: v.ppos ? dec(v.ppos,Float32Array) : null,
  scores:Object.fromEntries(Object.entries(v.scores).map(([k,b])=>[k,dec(b,Float32Array)])),
  S:v.structures.map(s=>({...s, nodes:dec(s.nodes,Float32Array), edges:dec(s.edges,Int32Array), emb:dec(s.emb,Float32Array), arc:dec(s.arc,Float32Array), curv:dec(s.curv,Float32Array), zeta:dec(s.zeta,Float32Array), degen:dec(s.degen,Float32Array), noise:dec(s.noise,Float32Array), cov:dec(s.cov,Float32Array), on:s.kind!=="group", colour:colourOf(s.id)}))}));
let view = "3d";
const TRI = D*(D+1)/2;
let lastThree = null;
const EXTRA = (DATA.extra||[]).map(e=>({name:e.name, kind:e.kind, cats:e.cats||[], values: e.kind==="num" ? dec(e.values,Float32Array) : dec(e.values,Int32Array)}));
function extraOf(key){ return key && key.startsWith("x:") ? EXTRA[+key.slice(2)] : null; }
function refreshExtraOptions(){
  const keep=(id,items)=>{ const el=$(id), old=el.value; fill(el, items); if(items.some(([v])=>v===old)) el.value=old; };
  const xs=EXTRA.map((e,i)=>["x:"+i, e.name+(e.kind==="cat"?" (categories)":"")]), xn=xs.filter((_,i)=>EXTRA[i].kind==="num");
  keep("colour", [["structure","structure"],["dim","dimension of the structure"]].concat(DATA.scores.map(s=>[s.key,s.label]), xs));
  keep("gcolour", [["count","points per node"],["median","median of the outlier score"],["nout","outliers per node"],["deg","node degree"],["curv","curvature at the node"],["noise","noise scale of the node"],["zeta","covariance scale zeta"],["degen","degenerate components"]].concat(DATA.names.map((n,i)=>["pc:"+i,"mean "+n]), xn.map(([v,t])=>[v,"mean "+t])));
  keep("gy", [["score","the outlier score"]].concat(DATA.names.map((n,i)=>["pc:"+i,n]), xn));
}
function pointValue(key){ // per-point numeric accessor for "pc:k" or a numeric extra
  if(key.startsWith("pc:")){ const k=+key.slice(3); return i=>X[i*D+k]; }
  const e=extraOf(key); return e && e.kind==="num" ? (i=>e.values[i]) : null; }
function parseLine(line, delim){ const out=[]; let cur="", q=false;
  for(let i=0;i<line.length;i++){ const c=line[i]; if(q){ if(c==='"'){ if(line[i+1]==='"'){cur+='"';i++;} else q=false; } else cur+=c; } else if(c==='"') q=true; else if(c===delim){ out.push(cur); cur=""; } else cur+=c; }
  out.push(cur); return out.map(t=>t.trim()); }
function loadExtraText(text, fname){
  const lines=text.split(/\r?\n/).filter(l=>l.trim().length); if(lines.length<2) throw new Error("the file has no data lines");
  const delim = lines[0].includes("\t") ? "\t" : lines[0].includes(";") && !lines[0].includes(",") ? ";" : ",";
  const head=parseLine(lines[0],delim), rows=lines.slice(1).map(l=>parseLine(l,delim));
  const ri=head.findIndex(h=>h.toLowerCase()==="row"); let index;
  if(ri>=0){ index=rows.map(r=>parseInt(r[ri],10)); if(index.some(k=>!(k>=0&&k<N))) throw new Error("the row column holds values outside 0.."+(N-1)); }
  else { if(rows.length!==N) throw new Error(`the file has ${rows.length} data lines but there are ${N} points; add a "row" column with the point index`); index=rows.map((_,k)=>k); }
  const added=[], skipped=[];
  head.forEach((name,c)=>{ if(c===ri||!name) return;
    const raw=new Array(N).fill(""); rows.forEach((r,k)=>{ raw[index[k]] = r[c]===undefined ? "" : r[c]; });
    const filled=raw.filter(t=>t!==""&&t.toLowerCase()!=="nan"); if(!filled.length){ skipped.push(name+" (empty)"); return; }
    let e;
    if(filled.every(t=>Number.isFinite(Number(t)))){ const vals=new Float32Array(N); for(let i=0;i<N;i++){ const t=raw[i]; vals[i] = (t===""||t.toLowerCase()==="nan") ? NaN : Number(t); } e={name,kind:"num",cats:[],values:vals}; }
    else { const cats=[...new Set(filled)].sort(); if(cats.length>60){ skipped.push(name+" (more than 60 categories)"); return; }
      const code=new Map(cats.map((t,k)=>[t,k])), vals=new Int32Array(N); for(let i=0;i<N;i++) vals[i]=code.has(raw[i])?code.get(raw[i]):-1; e={name,kind:"cat",cats,values:vals}; }
    const old=EXTRA.findIndex(x=>x.name===name); if(old>=0) EXTRA[old]=e; else EXTRA.push(e); added.push(name); });
  refreshExtraOptions();
  $("xinfo").textContent = `${fname}: ${rows.length.toLocaleString()} lines, added ${added.length? added.join(", "):"nothing"}`+(skipped.length? "; skipped "+skipped.join(", "):"");
  if(added.length){ const last=EXTRA.findIndex(x=>x.name===added[0]); $("colour").value="x:"+last; }
  draw();
}
const $ = id => document.getElementById(id);
let cur = 0, lastOut = null;
function fill(sel, items, val){ sel.innerHTML=""; items.forEach(([v,t])=>{const o=document.createElement("option");o.value=v;o.textContent=t;sel.appendChild(o);}); if(val!==undefined) sel.value=val; }
fill($("version"), V.map((v,i)=>[i,v.name]));
for(let a=0;a<3;a++) fill($("ax"+a), DATA.names.map((n,i)=>[i,n]).concat(a===2?[[-1,"(none: 2-D view)"]]:[]), a<D ? a : (a===2?-1:D-1));

fill($("score"), DATA.scores.map(s=>[s.key,s.label]));
function children(v, id){ return v.S.filter(s=>s.parent===id); }
function setSubtree(v, s, on){ s.on = on; children(v,s.id).forEach(c=>setSubtree(v,c,on)); }
function buildTree(){
  const v=V[cur], box=$("tree"); box.innerHTML=""; syncParents(v);
  const add=(s,depth)=>{
    const row=document.createElement("div"); row.className="node"; row.style.paddingLeft=(6+depth*16)+"px";
    const cb=document.createElement("input"); cb.type="checkbox"; cb.checked=s.on; cb.dataset.id=s.id;
    if(s.kind==="group"){ const st=leafStates(v,s); cb.indeterminate = st.some(x=>x) && !st.every(x=>x); }
    cb.addEventListener("change",()=>{ setSubtree(v,s,cb.checked); syncParents(v); buildTree(); draw(); });
    const sw=document.createElement("span"); sw.className="sw"; sw.style.background = s.kind==="group" ? "transparent" : s.colour; if(s.kind==="group") sw.style.border="1px solid var(--muted)";
    const lab=document.createElement("span"); lab.textContent=s.label;
    const kd=document.createElement("span"); kd.className="kind"; kd.textContent = s.dim+"-D"+(s.kind==="group"?" · parent":"");
    const cnt=document.createElement("span"); cnt.className="cnt"; cnt.textContent = s.kind==="group" ? subtreePoints(v,s).toLocaleString() : s.points.toLocaleString();
    row.append(cb,sw,lab,kd,cnt); box.appendChild(row);
    children(v,s.id).forEach(c=>add(c,depth+1));
  };
  children(v,-1).forEach(s=>add(s,0));
  $("summary").textContent=v.summary;
}
function leafStates(v,s){ return s.kind==="group" ? children(v,s.id).flatMap(c=>leafStates(v,c)) : [s.on]; }
function subtreePoints(v,s){ return s.kind==="group" ? children(v,s.id).reduce((a,c)=>a+subtreePoints(v,c),0) : s.points; }
function syncParents(v){ const rec=s=>{ if(s.kind!=="group") return s.on; const k=children(v,s.id).map(rec); s.on=k.some(x=>x); return s.on; }; children(v,-1).forEach(rec); }
function quantile(arr, q){ if(!arr.length) return NaN; const a=Float64Array.from(arr).sort(); const p=(a.length-1)*q, lo=Math.floor(p), hi=Math.ceil(p); return a[lo]+(a[hi]-a[lo])*(p-lo); }
function scoreMeta(){ return DATA.scores.find(s=>s.key===$("score").value); }
function setThrRange(){
  const v=V[cur], s=v.scores[$("score").value]; let lo=Infinity, hi=-Infinity;
  for(let i=0;i<N;i++){ const x=s[i]; if(Number.isFinite(x)){ if(x<lo)lo=x; if(x>hi)hi=x; } }
  if(!Number.isFinite(lo)){lo=0;hi=1;}
  const vals=[]; for(let i=0;i<N;i+=7){ if(Number.isFinite(s[i])) vals.push(s[i]); }
  const start = scoreMeta().low ? quantile(vals,0.01) : quantile(vals,0.99);
  $("thr").min=lo; $("thr").max=hi; $("thr").value=start; $("thrnum").value=(+start).toPrecision(4);
}
function compute(only){
  // only: if given, the selection is restricted to that one structure (graph view)
  const v=V[cur], on=v.S.map(s=> only===undefined ? (s.on && s.kind!=="group") : s.id===only), sid=v.sid, sc=v.scores[$("score").value], low=scoreMeta().low;
  const mode=document.querySelector("input[name=mode]:checked").value, scope=$("scope").value;
  const shown=new Uint8Array(N); let nshown=0;
  for(let i=0;i<N;i++){ const s=sid[i]; if(s>=0 && on[s]){shown[i]=1;nshown++;} }
  const out=new Uint8Array(N); let nout=0, info=""; const cuts=new Map();
  if(mode==="pct"){
    const q=Math.min(Math.max(+$("pctnum").value,0),100)/100;
    const groups = new Map();
    for(let i=0;i<N;i++){ if(!shown[i]||!Number.isFinite(sc[i])) continue; const key = scope==="each" ? sid[i] : -1; if(!groups.has(key)) groups.set(key,[]); groups.get(key).push(i); }
    const thr=[];
    groups.forEach((idx,key)=>{ const t=quantile(idx.map(i=>sc[i]), low?q:1-q); thr.push(t); cuts.set(key,t); idx.forEach(i=>{ if(low? sc[i]<t : sc[i]>t){out[i]=1;nout++;} }); });
    info = scope==="each" ? `${groups.size} structures, each with its own cut` : `cut at ${thr.length?thr[0].toPrecision(4):"–"}`;
  } else {
    const t=+$("thrnum").value; cuts.set(-1,t);
    for(let i=0;i<N;i++){ if(shown[i] && Number.isFinite(sc[i]) && (low? sc[i]<t : sc[i]>t)){out[i]=1;nout++;} }
    info=`${low?"score below":"score above"} ${t.toPrecision(4)}`;
  }
  return {shown,out,nshown,nout,info,low,cuts};
}
function seqColour(vals){ // single-hue ramp, robust range
  const f=vals.filter(Number.isFinite); const lo=quantile(f,0.02), hi=quantile(f,0.98); return {cmin:lo,cmax:hi}; }
function draw(){ buildDimToggles(); catLegend(); if(view==="graph") drawGraph(); else draw3d(); }
function catLegend(){ const e=extraOf($("colour").value), box=$("catlegend"); box.innerHTML="";
  if(!e||e.kind!=="cat"||view==="graph") return;
  e.cats.forEach((t,k)=>{ const sp=document.createElement("span"), sw=document.createElement("span"); sw.className="sw"; sw.style.background=colourOf(k); sp.append(sw, document.createTextNode(t)); box.appendChild(sp); }); }
function buildDimToggles(){
  const v=V[cur], box=$("dimtog"); box.innerHTML="";
  const dims=[...new Set(v.S.filter(s=>s.kind!=="group").map(s=>s.dim))].sort((a,b)=>a-b);
  dims.forEach(dm=>{
    const mine=v.S.filter(s=>s.kind!=="group" && s.dim===dm), st=mine.map(s=>s.on);
    const lab=document.createElement("label"), cb=document.createElement("input"); cb.type="checkbox";
    cb.checked=st.every(x=>x); cb.indeterminate=st.some(x=>x)&&!st.every(x=>x);
    cb.addEventListener("change",()=>{ mine.forEach(s=>{s.on=cb.checked;}); syncParents(v); buildTree(); draw(); });
    const n=mine.reduce((a,s)=>a+s.points,0);
    lab.append(cb, document.createTextNode(`${dm}-D`), Object.assign(document.createElement("span"),{className:"kind",textContent:`(${mine.length}, ${n.toLocaleString()} pts)`}));
    box.appendChild(lab); });
}
function fillGraphSelect(){
  const v=V[cur], opts=v.S.filter(s=>s.kind!=="group" && s.points>0).sort((a,b)=>b.points-a.points).map(s=>[s.id, `${s.label} · ${s.dim}-D · ${s.points.toLocaleString()} points`]);
  fill($("gstruct"), opts);
}
function hash01(i,k){ let x=((i+1)*2654435761 + k*40503)>>>0; x^=x>>>15; x=Math.imul(x,2246822519)>>>0; x^=x>>>13; return x/4294967296; }
function median(a){ if(!a.length) return NaN; const b=Float64Array.from(a).sort(); const m=b.length>>1; return b.length%2? b[m] : 0.5*(b[m-1]+b[m]); }
const RAMP=["#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#08519c","#08306b"];
const RAMPSCALE=RAMP.map((c,i)=>[i/(RAMP.length-1),c]);
function ramp(vals){ const f=Array.from(vals).filter(Number.isFinite); if(!f.length) return {lo:0,hi:1,bin:()=>-1};
  let lo=quantile(f,0.02), hi=quantile(f,0.98); if(!(hi>lo)){ lo=Math.min(...f); hi=Math.max(...f); if(!(hi>lo)) hi=lo+1e-9; }
  return {lo,hi,bin:x=>Number.isFinite(x)? Math.max(0,Math.min(RAMP.length-1,Math.floor((x-lo)/(hi-lo)*RAMP.length))) : -1}; }
function nodeProperty(v,g,idx,r,sc,cm,meta){
  const K=g.n_nodes, E=g.edges, n=E.length/2, cnt=new Float64Array(K), nout=new Float64Array(K), val=new Array(K).fill(NaN); let title="";
  idx.forEach(i=>{ cnt[v.node[i]]++; if(r.out[i]) nout[v.node[i]]++; });
  if(cm==="count"){ for(let k=0;k<K;k++) val[k]=cnt[k]; title="points per node"; }
  else if(cm==="nout"){ for(let k=0;k<K;k++) val[k]=nout[k]; title="outliers per node"; }
  else if(cm==="median"){ const lists=Array.from({length:K},()=>[]); idx.forEach(i=>{ if(Number.isFinite(sc[i])) lists[v.node[i]].push(sc[i]); }); for(let k=0;k<K;k++) val[k]=median(lists[k]); title="median "+meta.label; }
  else if(cm==="deg"){ val.fill(0); for(let e=0;e<n;e++){ val[E[2*e]]++; val[E[2*e+1]]++; } title="node degree (number of edges)"; }
  else if(cm==="curv"){ const sum=new Float64Array(K), c=new Float64Array(K); for(let e=0;e<n;e++){ const q=g.curv[e]; if(!Number.isFinite(q)) continue; sum[E[2*e]]+=q; c[E[2*e]]++; sum[E[2*e+1]]+=q; c[E[2*e+1]]++; } for(let k=0;k<K;k++) val[k]=c[k]? sum[k]/c[k] : NaN; title="mean curvature of the edges at the node"; }
  else if(cm==="noise"){ for(let k=0;k<K;k++) val[k]=g.noise[k]; title="noise scale of the node's Gaussian"; }
  else if(cm==="zeta"){ for(let k=0;k<K;k++) val[k]=g.zeta[k]; title="covariance scale zeta"; }
  else if(cm==="degen"){ for(let k=0;k<K;k++) val[k]=g.degen[k]; title="degenerate component (1 = yes)"; }
  else { const pf=pointValue(cm), sum=new Float64Array(K), c=new Float64Array(K); idx.forEach(i=>{ const t=pf(i); if(Number.isFinite(t)){ sum[v.node[i]]+=t; c[v.node[i]]++; } });
    for(let k=0;k<K;k++) val[k]=c[k]? sum[k]/c[k] : NaN; title="mean "+(cm.startsWith("pc:")? DATA.names[+cm.slice(3)] : extraOf(cm).name); }
  return {val,title,cnt,nout};
}
function edgeProperty(g,key,pos,nd){
  const E=g.edges, n=E.length/2, P=g.nodes, val=new Float64Array(n); let title="";
  const elen=e=>{ let q=0; const a=E[2*e]*D, b=E[2*e+1]*D; for(let k=0;k<D;k++){ const t=P[a+k]-P[b+k]; q+=t*t; } return Math.sqrt(q); };
  if(key==="curv"){ for(let e=0;e<n;e++) val[e]=g.curv[e]; title="curvature along the edge"; }
  else if(key==="len"){ for(let e=0;e<n;e++) val[e]=elen(e); title="edge length in data space"; }
  else if(key==="stretch"){ const ratio=[]; for(let e=0;e<n;e++){ let q=0; for(let a=0;a<nd;a++){ const t=pos(E[2*e],a)-pos(E[2*e+1],a); q+=t*t; } val[e]=Math.sqrt(q)/Math.max(elen(e),1e-300); ratio.push(val[e]); }
    const md=median(ratio)||1; for(let e=0;e<n;e++) val[e]/=md; title="distortion of the layout (drawn length / true length, 1 = typical)"; }
  return {val,title};
}
function drawGraph(){
  const v=V[cur], sidSel=+$("gstruct").value, s=v.S[sidSel];
  if(!s){ Plotly.react("plot",[],{}); return; }
  const g=v.S[s.graph_of], K=g.n_nodes, E=g.edges, n=E.length/2, P=g.nodes, m=g.m;
  const r=compute(sidSel); lastOut=r;
  const key=$("score").value, sc=v.scores[key], meta=scoreMeta(), hl=$("hl").checked, only=$("only").checked, showPts=$("gpts").checked, bySize=$("gsize").checked;
  const idx=[]; for(let i=0;i<N;i++) if(v.sid[i]===sidSel) idx.push(i);
  const ht="row %{customdata[0]}<br>score %{customdata[1]:.4g}<extra></extra>";
  const nht="node %{customdata[0]}<br>%{customdata[1]} points, %{customdata[2]} outliers<br>value %{customdata[3]:.4g}<extra></extra>";
  const is1d = g.dim===1 && g.arc.length===K;
  const cutKey = (document.querySelector("input[name=mode]:checked").value==="pct" && $("scope").value==="each") ? sidSel : -1, cut = r.cuts.has(cutKey) ? r.cuts.get(cutKey) : NaN;
  const np=nodeProperty(v,g,idx,r,sc,$("gcolour").value,meta), nr=ramp(np.val);
  const maxc=Math.max(1,...np.cnt), nsize=k=> bySize ? 5+13*Math.sqrt(np.cnt[k]/maxc) : 9;
  const ekey=$("gedge").value;
  $("gprojwrap").style.display = (!is1d && g.dim>=4 && m>=4) ? "block" : "none";
  const traces=[]; let layout, kind;
  const nodeMarker=(sizes)=>({size:sizes,color:np.val,colorscale:RAMPSCALE,cmin:nr.lo,cmax:nr.hi,line:{color:"#ffffff",width:0.6},
    colorbar:{title:{text:np.title,side:"right"},thickness:12,len:0.55,x:1.0,y:0.72}});
  const edgeBar=(ep,er,type,extra)=>Object.assign({type,mode:"markers",x:[null],y:[null],hoverinfo:"skip",showlegend:false,
    marker:{size:0.1,color:[er.lo,er.hi],colorscale:RAMPSCALE,cmin:er.lo,cmax:er.hi,showscale:true,colorbar:{title:{text:ep.title,side:"right"},thickness:12,len:0.4,x:1.0,y:0.2}}},extra||{});
  if(is1d){
    kind="curve";
    const pos=(k,a)=> a===0 ? g.arc[k] : 0, ep = ekey==="none" ? null : edgeProperty(g,ekey,pos,1), er = ep ? ramp(ep.val) : null;
    const elen=(a,b)=>{ let q=0; for(let k=0;k<D;k++){ const t=P[a*D+k]-P[b*D+k]; q+=t*t; } return Math.sqrt(q); };
    const px=[],py=[],cd=[],ox=[],oy=[],ocd=[];
    const yk=$("gy").value, yf = yk==="score" ? (i=>sc[i]) : pointValue(yk), ylab = yk==="score" ? meta.label : (yk.startsWith("pc:") ? DATA.names[+yk.slice(3)] : extraOf(yk).name);
    idx.forEach(i=>{ const a=v.node[i], b=v.toward[i]; const x=g.arc[a]+Math.sign(g.arc[b]-g.arc[a])*v.frac[i]*elen(a,b), y=yf(i);
      if(!Number.isFinite(y)) return;
      if(r.out[i]){ ox.push(x);oy.push(y);ocd.push([i,sc[i]]); if(hl||only) return; } if(only) return; px.push(x);py.push(y);cd.push([i,sc[i]]); });
    if(showPts) traces.push({type:"scattergl",mode:"markers",x:px,y:py,customdata:cd,hovertemplate:ht,marker:{size:5,color:s.colour,opacity:0.7},showlegend:false});
    if(ox.length) traces.push({type:"scattergl",mode:"markers",x:ox,y:oy,customdata:ocd,hovertemplate:ht,marker:{size:9,color:"#e34948",symbol:"diamond",line:{color:"#fff",width:0.5}},showlegend:false});
    // the curve itself as a line: one segment per edge, nodes on top
    for(let e=0;e<n;e++){ const b=ep? er.bin(ep.val[e]) : -1;
      traces.push({type:"scatter",mode:"lines",x:[g.arc[E[2*e]],g.arc[E[2*e+1]]],y:[0,0],xaxis:"x",yaxis:"y2",line:{color:b>=0?RAMP[b]:"#b9b8b2",width:6},hoverinfo:"skip",showlegend:false}); }
    if(ep) traces.push(edgeBar(ep,er,"scatter",{xaxis:"x",yaxis:"y2"}));
    const nx=[],ny=[],ncd=[],sz=[]; for(let k=0;k<K;k++){ nx.push(g.arc[k]); ny.push(0); ncd.push([k,np.cnt[k],np.nout[k],np.val[k]]); sz.push(nsize(k)+3); }
    traces.push({type:"scatter",mode:"markers",x:nx,y:ny,xaxis:"x",yaxis:"y2",customdata:ncd,hovertemplate:nht,marker:nodeMarker(sz),showlegend:false});
    const shapes=[]; for(let k=0;k<K;k++) shapes.push({type:"line",xref:"x",yref:"y domain",x0:g.arc[k],x1:g.arc[k],y0:0,y1:1,line:{color:"#eeede8",width:1}});
    if(Number.isFinite(cut) && yk==="score") shapes.push({type:"line",xref:"x domain",yref:"y",x0:0,x1:1,y0:cut,y1:cut,line:{color:"#e34948",width:1,dash:"dash"}});
    layout={shapes, xaxis:{title:{text:"position along the curve (arc length)"},gridcolor:"#eeede8",zeroline:false,anchor:"y2"},
      yaxis:{title:{text:ylab},gridcolor:"#eeede8",zeroline:false,domain:[0.3,1]}, yaxis2:{domain:[0,0.16],range:[-1,1],visible:false,fixedrange:true},
      annotations:[{text:"the curve: nodes and edges",xref:"paper",yref:"paper",x:0,y:0.2,showarrow:false,font:{size:11,color:"#52514e"},xanchor:"left"}],
      margin:{l:60,r:90,t:10,b:50},dragmode:"pan"};
  } else {
    const three = g.dim>=3 && m>=3; kind = three ? (g.dim>=4? "projection of a "+g.dim+"-D structure" : "3-D layout") : "planar layout";
    let ax=[0,1,2]; if(three && g.dim>=4 && m>=4){ ax=$("gproj").value.split("").map(t=>+t-1); }
    const nd=three?3:2, pos=(k,a)=>g.emb[k*m+ax[a]], T=three?"scatter3d":"scattergl";
    const ep = ekey==="none" ? null : edgeProperty(g,ekey,pos,nd), er = ep ? ramp(ep.val) : null;
    const groups = ep ? RAMP.map(()=>[]) : [[]]; const none=[];
    for(let e=0;e<n;e++){ const b=ep? er.bin(ep.val[e]) : 0; (b>=0? groups[b] : none).push(e); }
    const lineTrace=(list,color,width)=>{ const xs=[],ys=[],zs=[]; list.forEach(e=>{ const a=E[2*e], b=E[2*e+1]; xs.push(pos(a,0),pos(b,0),null); ys.push(pos(a,1),pos(b,1),null); if(three) zs.push(pos(a,2),pos(b,2),null); });
      const t={type:T,mode:"lines",x:xs,y:ys,line:{color,width},hoverinfo:"skip",showlegend:false}; if(three) t.z=zs; return t; };
    if(none.length) traces.push(lineTrace(none,"#dddcd6",three?2:1));
    groups.forEach((list,b)=>{ if(list.length) traces.push(lineTrace(list, ep? RAMP[b] : "#b9b8b2", ep? (three?4:2) : (three?2:1))); });
    if(ep){ const eb=edgeBar(ep,er,T); if(three) eb.z=[null]; traces.push(eb); }
    // points at their position on the graph
    const med=[]; for(let e=0;e<Math.min(n,2000);e++){ let q=0; for(let a=0;a<nd;a++){ const t=pos(E[2*e],a)-pos(E[2*e+1],a); q+=t*t; } med.push(Math.sqrt(q)); }
    const jit=0.25*(median(med)||0.01), fill=$("gfill").checked && v.ppos;
    const px=[],py=[],pz=[],cd=[],ox=[],oy=[],oz=[],ocd=[];
    idx.forEach(i=>{ const a=v.node[i], b=v.toward[i], f=v.frac[i], c=[];
      for(let q=0;q<nd;q++) c.push(fill ? v.ppos[i*4+ax[q]] : pos(a,q)+f*(pos(b,q)-pos(a,q))+jit*(hash01(i,q+1)-0.5)*1.6);
      if(r.out[i]){ ox.push(c[0]);oy.push(c[1]);if(three)oz.push(c[2]);ocd.push([i,sc[i]]); if(hl||only) return; } if(only) return; px.push(c[0]);py.push(c[1]);if(three)pz.push(c[2]);cd.push([i,sc[i]]); });
    if(showPts && px.length){ const t={type:T,mode:"markers",x:px,y:py,customdata:cd,hovertemplate:ht,marker:{size:three?1.6:3,color:"#b9b8b2",opacity:three?0.5:0.6,line:{width:0}},showlegend:false}; if(three) t.z=pz; traces.push(t); }
    const nx=[],ny=[],nz=[],ncd=[],sz=[]; for(let k=0;k<K;k++){ nx.push(pos(k,0)); ny.push(pos(k,1)); if(three) nz.push(pos(k,2)); ncd.push([k,np.cnt[k],np.nout[k],np.val[k]]); sz.push(three? nsize(k)*0.55 : nsize(k)); }
    const nt={type:T,mode:"markers",x:nx,y:ny,customdata:ncd,hovertemplate:nht,marker:nodeMarker(sz),showlegend:false}; if(three){ nt.z=nz; nt.marker.line={width:0}; } traces.push(nt);
    if(ox.length){ const t={type:T,mode:"markers",x:ox,y:oy,customdata:ocd,hovertemplate:ht,marker:{size:three?4.5:9,color:"#e34948",symbol:"diamond",line:{color:"#fff",width:0.5}},showlegend:false}; if(three) t.z=oz; traces.push(t); }
    if(three){ const axn=a=>({title:{text:"intrinsic coordinate "+(ax[a]+1)},backgroundcolor:"#ffffff",gridcolor:"#eeede8",zerolinecolor:"#eeede8",color:"#52514e",showticklabels:false});
      layout={scene:{aspectmode:"data",xaxis:axn(0),yaxis:axn(1),zaxis:axn(2)},margin:{l:0,r:0,t:0,b:0}}; }
    else layout={xaxis:{visible:false},yaxis:{visible:false,scaleanchor:"x"},margin:{l:10,r:90,t:10,b:10},dragmode:"pan"};
  }
  Plotly.react("plot", traces, Object.assign({paper_bgcolor:"#fcfcfb",plot_bgcolor:"#ffffff",uirevision:"g"+cur+"-"+sidSel+"-"+$("gproj").value},layout), {displaylogo:false,responsive:true,scrollZoom:true});
  drawLocator(v,s,idx,r);
  $("outinfo").innerHTML = `<b>${r.nout.toLocaleString()}</b> outliers among ${r.nshown.toLocaleString()} points of this structure (${r.nshown? (100*r.nout/r.nshown).toFixed(2):"0"} %)<br><span class="kind">${meta.low?"low":"high"} values are outliers · ${Number.isFinite(cut)? "cut at "+cut.toPrecision(4) : r.info}</span>`;
  $("status").textContent = `${v.name} · structure ${s.label}: ${s.dim}-D, ${K} nodes, ${n} edges, ${idx.length.toLocaleString()} points · ${kind}`;
}
function drawLocator(v,s,idx,r){
  // inset: the selected structure among all the data, on the features chosen under "View"
  const box=$("loc"); box.classList.toggle("off", !$("gloc").checked); if(!$("gloc").checked) return;
  const three=$("ax2").value!=="-1", ax = three ? [+$("ax0").value,+$("ax1").value,+$("ax2").value] : [+$("ax0").value,+$("ax1").value], nd=ax.length, T=three?"scatter3d":"scattergl";
  if(box.dataset.three!==String(three)){ Plotly.purge("locplot"); box.dataset.three=String(three); }
  const mk=()=>ax.map(()=>[]), B=mk(), S=mk(), O=mk(), bstep=Math.max(1,Math.ceil(N/8000)), sstep=Math.max(1,Math.ceil(idx.length/6000));
  for(let i=0;i<N;i+=bstep){ if(v.sid[i]===s.id) continue; for(let a=0;a<nd;a++) B[a].push(X[i*D+ax[a]]); }
  idx.forEach((i,q)=>{ if(r.out[i]){ for(let a=0;a<nd;a++) O[a].push(X[i*D+ax[a]]); } else if(q%sstep===0){ for(let a=0;a<nd;a++) S[a].push(X[i*D+ax[a]]); } });
  const tr=(C,marker)=>{ const t={type:T,mode:"markers",x:C[0],y:C[1],marker,hoverinfo:"skip",showlegend:false}; if(three) t.z=C[2]; return t; };
  const traces=[tr(B,{size:three?1.2:2,color:"#c9c8c2",opacity:three?0.35:0.4,line:{width:0}}), tr(S,{size:three?2:3,color:s.colour,opacity:0.9,line:{width:0}})];
  if(O[0].length) traces.push(tr(O,{size:three?2.2:5,color:"#e34948",symbol:"diamond",line:{width:0}}));
  const axn=a=>({title:{text:DATA.names[ax[a]],font:{size:10}},tickfont:{size:8},backgroundcolor:"#ffffff",gridcolor:"#eeede8",zerolinecolor:"#eeede8",color:"#52514e"});
  const layout = three ? {scene:{aspectmode:"data",xaxis:axn(0),yaxis:axn(1),zaxis:axn(2)},margin:{l:0,r:0,t:0,b:0}}
    : {xaxis:Object.assign(axn(0),{zeroline:false}),yaxis:Object.assign(axn(1),{zeroline:false,scaleanchor:"x"}),margin:{l:38,r:6,t:6,b:32},dragmode:"pan"};
  Plotly.react("locplot",traces,Object.assign({paper_bgcolor:"#ffffff",plot_bgcolor:"#ffffff",uirevision:"loc"+three},layout),{displayModeBar:false,responsive:true,scrollZoom:true});
  $("loctitle").textContent = `structure ${s.label} in the data (${ax.map(a=>DATA.names[a]).join(", ")})`;
}
function setView(w){ view=w; document.body.className = w==="graph" ? "vgraph" : "v3d"; document.querySelectorAll(".tab").forEach(t=>t.classList.toggle("active", t.dataset.view===w)); Plotly.purge("plot"); lastThree=null; draw(); }
function triIndex(a,b){ if(a>b){ const t=a;a=b;b=t; } return a*D - (a*(a-1))/2 + (b-a); }
function gridDensity(structs, ax, G){
  // marginal density of the selected leaves on the shown features, splatted on a regular grid
  const nd=ax.length, lo=new Array(nd).fill(Infinity), hi=new Array(nd).fill(-Infinity), v=V[cur];
  structs.forEach(s=>{ const g=v.S[s.graph_of], K=g.n_nodes; for(let k=0;k<K;k++) for(let a=0;a<nd;a++){ const mu=g.nodes[k*D+ax[a]], sd=Math.sqrt(Math.max(g.cov[k*TRI+triIndex(ax[a],ax[a])],0)); if(mu-3.5*sd<lo[a]) lo[a]=mu-3.5*sd; if(mu+3.5*sd>hi[a]) hi[a]=mu+3.5*sd; } });
  if(!lo.every(Number.isFinite)) return null;
  const h=lo.map((l,a)=>(hi[a]-l)/(G-1)), n=Math.pow(G,nd), val=new Float32Array(n); let wsum=0; structs.forEach(s=>wsum+=s.weight);
  structs.forEach(s=>{ const g=v.S[s.graph_of], K=g.n_nodes, w=s.weight/wsum/K;
    for(let k=0;k<K;k++){
      const mu=[], C=[]; for(let a=0;a<nd;a++){ mu.push(g.nodes[k*D+ax[a]]); for(let b=0;b<nd;b++) C.push(g.cov[k*TRI+triIndex(ax[a],ax[b])]); }
      let inv, det;
      if(nd===2){ det=C[0]*C[3]-C[1]*C[2]; if(!(det>0)) continue; inv=[C[3]/det,-C[1]/det,-C[2]/det,C[0]/det]; }
      else { const c00=C[4]*C[8]-C[5]*C[7], c01=C[5]*C[6]-C[3]*C[8], c02=C[3]*C[7]-C[4]*C[6]; det=C[0]*c00+C[1]*c01+C[2]*c02; if(!(det>0)) continue;
        inv=[c00/det,(C[2]*C[7]-C[1]*C[8])/det,(C[1]*C[5]-C[2]*C[4])/det, c01/det,(C[0]*C[8]-C[2]*C[6])/det,(C[2]*C[3]-C[0]*C[5])/det, c02/det,(C[1]*C[6]-C[0]*C[7])/det,(C[0]*C[4]-C[1]*C[3])/det]; }
      const norm=w/Math.sqrt(Math.pow(2*Math.PI,nd)*det), i0=[], i1=[];
      for(let a=0;a<nd;a++){ const sd=Math.sqrt(C[a*nd+a]); i0.push(Math.max(0,Math.floor((mu[a]-3.5*sd-lo[a])/h[a]))); i1.push(Math.min(G-1,Math.ceil((mu[a]+3.5*sd-lo[a])/h[a]))); }
      if(nd===2){ for(let j=i0[1];j<=i1[1];j++){ const dy=lo[1]+j*h[1]-mu[1]; for(let i=i0[0];i<=i1[0];i++){ const dx=lo[0]+i*h[0]-mu[0]; const q=dx*(inv[0]*dx+inv[1]*dy)+dy*(inv[2]*dx+inv[3]*dy); if(q<25) val[j*G+i]+=norm*Math.exp(-0.5*q); } } }
      else { for(let l=i0[2];l<=i1[2];l++){ const dz=lo[2]+l*h[2]-mu[2]; for(let j=i0[1];j<=i1[1];j++){ const dy=lo[1]+j*h[1]-mu[1]; for(let i=i0[0];i<=i1[0];i++){ const dx=lo[0]+i*h[0]-mu[0];
        const q=dx*(inv[0]*dx+inv[1]*dy+inv[2]*dz)+dy*(inv[3]*dx+inv[4]*dy+inv[5]*dz)+dz*(inv[6]*dx+inv[7]*dy+inv[8]*dz); if(q<25) val[(l*G+j)*G+i]+=norm*Math.exp(-0.5*q); } } } }
    } });
  return {lo,h,G,nd,val};
}
function massLevels(grid, masses){ // density value whose super-level set holds the given share of the model's mass
  const s=Float32Array.from(grid.val).sort(), n=s.length; let tot=0; for(let i=0;i<n;i++) tot+=s[i];
  const out=[]; masses.forEach(p=>{ let acc=0, lvl=s[n-1]; for(let i=n-1;i>=0;i--){ acc+=s[i]; if(acc>=p*tot){ lvl=s[i]; break; } } out.push(lvl); }); return out; }
function densityTraces(ax){
  const mode=$("isomode").value; if(mode==="off") return {traces:[],note:""};
  const v=V[cur], three=ax.length===3, masses=[...document.querySelectorAll(".isolvl:checked")].map(e=>+e.value).sort((a,b)=>a-b);
  if(!masses.length) return {traces:[],note:""};
  let sel=v.S.filter(s=>s.on && s.kind!=="group" && v.S[s.graph_of].n_nodes>0), note="";
  if(!sel.length) return {traces:[],note:""};
  const G = three ? +$("isores").value : 4*(+$("isores").value);
  const groups = mode==="each" ? sel.sort((a,b)=>b.weight-a.weight).slice(0,12).map(s=>[s]) : [sel];
  if(mode==="each" && sel.length>12) note=" · density drawn for the 12 heaviest ticked structures";
  const traces=[];
  groups.forEach(gr=>{ const grid=gridDensity(gr,ax,G); if(!grid) return; const lv=massLevels(grid,masses), col = mode==="each" ? gr[0].colour : "#17457f";
    const coords=a=>Array.from({length:G},(_,i)=>grid.lo[a]+i*grid.h[a]);
    if(three){ const n=G*G*G, xs=new Float32Array(n), ys=new Float32Array(n), zs=new Float32Array(n); let t=0;
      for(let l=0;l<G;l++) for(let j=0;j<G;j++) for(let i=0;i<G;i++){ xs[t]=grid.lo[0]+i*grid.h[0]; ys[t]=grid.lo[1]+j*grid.h[1]; zs[t]=grid.lo[2]+l*grid.h[2]; t++; }
      lv.forEach((L,q)=>{ if(!(L>0)) return; traces.push({type:"isosurface",x:xs,y:ys,z:zs,value:grid.val,isomin:L,isomax:L,surface:{count:1},caps:{x:{show:false},y:{show:false},z:{show:false}},
        colorscale:[[0,col],[1,col]],showscale:false,opacity: masses[q]<=0.5 ? 0.35 : masses[q]<=0.9 ? 0.22 : masses[q]<=0.99 ? 0.13 : 0.08,hoverinfo:"skip",lighting:{ambient:0.8,diffuse:0.5,specular:0.1},name:"density"}); });
    } else { const z=[]; for(let j=0;j<G;j++) z.push(Array.from(grid.val.subarray(j*G,(j+1)*G)));
      lv.forEach((L,q)=>{ if(!(L>0)) return; traces.push({type:"contour",x:coords(0),y:coords(1),z,autocontour:false,contours:{coloring:"none",start:L,end:L,size:1,showlabels:false},
        line:{color:col,width:[2.2,1.6,1.1,0.8][Math.min(q,3)],dash:["solid","solid","dash","dot"][Math.min(q,3)]},showscale:false,showlegend:false,hoverinfo:"skip",name:"density"}); }); }
  });
  return {traces, note:note+` · density ${three?"surfaces":"curves"} enclose ${masses.map(m=>(100*m)+" %").join(", ")} of the model's mass in this projection`};
}
function draw3d(){
  const v=V[cur], cm=$("colour").value, three = $("ax2").value!=="-1";
  const ax = three ? [+$("ax0").value,+$("ax1").value,+$("ax2").value] : [+$("ax0").value,+$("ax1").value], nd=ax.length, T=three?"scatter3d":"scattergl";
  if(lastThree!==three){ Plotly.purge("plot"); lastThree=three; }
  const r=compute(); lastOut=r;
  const only=$("only").checked, hl=$("hl").checked, sc=v.scores[$("score").value], xe=extraOf(cm);
  const P=[[],[],[]], O=[[],[],[]], pc=[], cd=[], ocd=[];
  for(let i=0;i<N;i++){
    if(!r.shown[i]) continue; const b=i*D, isOut=r.out[i]===1;
    if(isOut && (hl||only)){ for(let a=0;a<nd;a++) O[a].push(X[b+ax[a]]); ocd.push([i,v.S[v.sid[i]].label,sc[i]]); if(hl) continue; }
    if(only) continue;
    for(let a=0;a<nd;a++) P[a].push(X[b+ax[a]]); cd.push([i,v.S[v.sid[i]].label,sc[i]]);
    if(cm==="structure") pc.push(v.S[v.sid[i]].colour); else if(cm==="dim") pc.push(DIMCOL[Math.min(v.S[v.sid[i]].dim,8)]);
    else if(xe) pc.push(xe.kind==="num" ? xe.values[i] : (xe.values[i]<0 ? "#b9b8b2" : colourOf(xe.values[i]))); else pc.push(v.scores[cm][i]);
  }
  const traces=[], put=(t,A)=>{ t.x=A[0]; t.y=A[1]; if(three) t.z=A[2]; return t; };
  const ht="row %{customdata[0]} · %{customdata[1]}<br>score %{customdata[2]:.4g}<extra></extra>";
  const dens=densityTraces(ax);
  if(!three) dens.traces.forEach(t=>traces.push(t));
  if($("showpts").checked && P[0].length){
    const marker={size:+$("psize").value+(three?0:1), opacity:+$("popac").value, color:pc, line:{width:0}};
    if(cm!=="structure" && cm!=="dim" && !(xe && xe.kind==="cat")){ const rg=seqColour(pc); Object.assign(marker,{colorscale:"Blues",reversescale:true,cmin:rg.cmin,cmax:rg.cmax,
      colorbar:{title:{text: xe ? xe.name : DATA.scores.find(s=>s.key===cm).label,side:"right"},thickness:12,len:0.6}}); }
    traces.push(put({type:T,mode:"markers",customdata:cd,hovertemplate:ht,marker,name:"points",showlegend:false},P));
  }
  if($("showgraphs").checked){
    const done=new Set();
    v.S.forEach(s=>{ if(!s.on||s.kind==="group") return; const g=v.S[s.graph_of]; if(done.has(g.id)||!g.n_nodes) return; done.add(g.id);
      const E=g.edges, Q=g.nodes, n=E.length/2, L=[new Array(3*n),new Array(3*n),new Array(3*n)];
      for(let e=0;e<n;e++){ const a=E[2*e]*D, b=E[2*e+1]*D; for(let q=0;q<nd;q++){ L[q][3*e]=Q[a+ax[q]]; L[q][3*e+1]=Q[b+ax[q]]; L[q][3*e+2]=null; } }
      traces.push(put({type:T,mode:"lines",line:{color: cm==="structure"? (s.kind==="residual"? s.colour : g.colour) : "#0b0b0b", width:three?3:1.2},hoverinfo:"skip",name:"graph "+g.label,showlegend:false},L)); });
  }
  if(three) dens.traces.forEach(t=>traces.push(t));
  if(O[0].length) traces.push(put({type:T,mode:"markers",customdata:ocd,hovertemplate:ht,marker:{size:+$("psize").value+(three?2.5:5),color:"#e34948",symbol:"diamond",line:{color:"#ffffff",width:0.5}},name:"outliers",showlegend:false},O));
  const axis=t=>({title:{text:t},backgroundcolor:"#ffffff",gridcolor:"#dddcd6",zerolinecolor:"#dddcd6",color:"#52514e"});
  const layout = three ? {margin:{l:0,r:0,t:0,b:0},paper_bgcolor:"#fcfcfb",uirevision:"keep3",scene:{aspectmode:"data",xaxis:axis(DATA.names[ax[0]]),yaxis:axis(DATA.names[ax[1]]),zaxis:axis(DATA.names[ax[2]]),uirevision:"keep3"}}
    : {margin:{l:60,r:20,t:10,b:50},paper_bgcolor:"#fcfcfb",plot_bgcolor:"#ffffff",uirevision:"keep2",dragmode:"pan",xaxis:{title:{text:DATA.names[ax[0]]},gridcolor:"#eeede8",zeroline:false},yaxis:{title:{text:DATA.names[ax[1]]},gridcolor:"#eeede8",zeroline:false,scaleanchor:"x"}};
  Plotly.react("plot", traces, layout, {displaylogo:false,responsive:true,scrollZoom:!three});
  $("outinfo").innerHTML = `<b>${r.nout.toLocaleString()}</b> outliers among ${r.nshown.toLocaleString()} shown points (${r.nshown? (100*r.nout/r.nshown).toFixed(2):"0"} %)<br><span class="kind">${scoreMeta().low?"low":"high"} values are outliers · ${r.info}</span>`;
  $("status").textContent = `${v.name} · ${r.nshown.toLocaleString()} of ${N.toLocaleString()} points shown${dens.note} · ${three?"drag to rotate, scroll to zoom":"drag to pan, scroll to zoom"}`;
}
function download(onlyOut){
  const v=V[cur], r=lastOut||compute(); const keys=DATA.scores.map(s=>s.key);
  const rows=["row,structure,dim,node,outlier,"+keys.join(",")+","+DATA.names.join(",")];
  for(let i=0;i<N;i++){ if(!r.shown[i]||(onlyOut&&!r.out[i])) continue; const s=v.S[v.sid[i]];
    const xs=[]; for(let k=0;k<D;k++) xs.push(X[i*D+k].toPrecision(6));
    rows.push([i,'"'+s.label+'"',s.dim,v.node[i],r.out[i],...keys.map(k=>v.scores[k][i]),...xs].join(",")); }
  const a=document.createElement("a"); a.href=URL.createObjectURL(new Blob([rows.join("\n")],{type:"text/csv"})); a.download=(onlyOut?"outliers":"points")+".csv"; a.click(); URL.revokeObjectURL(a.href);
}
function setAll(fn){ const v=V[cur]; v.S.forEach(s=>{ s.on=fn(s,v); }); syncParents(v); buildTree(); draw(); }
$("all").onclick=()=>setAll(()=>true); $("none").onclick=()=>setAll(()=>false);
$("leavesonly").onclick=()=>setAll(s=>s.kind==="leaf");
$("top").onclick=()=>setAll(s=>s.kind==="residual");
$("version").onchange=()=>{ cur=+$("version").value; buildTree(); fillGraphSelect(); setThrRange(); draw(); };
document.querySelectorAll(".tab").forEach(t=>t.addEventListener("click",()=>setView(t.dataset.view)));
refreshExtraOptions();
["gstruct","gcolour","gpts","gfill","gloc","gy","gedge","gproj","gsize","isomode","isores"].forEach(id=>$(id).addEventListener("input",draw));
document.querySelectorAll(".isolvl").forEach(e=>e.addEventListener("change",draw));
$("xfile").addEventListener("change",ev=>{ const f=ev.target.files[0]; if(!f) return; const rd=new FileReader();
  rd.onload=()=>{ try{ loadExtraText(String(rd.result), f.name); }catch(err){ $("xinfo").textContent="Could not use "+f.name+": "+err.message; } };
  rd.readAsText(f); });
["ax0","ax1","ax2","colour","psize","popac","showpts","showgraphs","scope","hl","only"].forEach(id=>$(id).addEventListener("input",draw));
$("score").onchange=()=>{ setThrRange(); draw(); };
$("pct").oninput=()=>{ $("pctnum").value=$("pct").value; document.querySelector("input[name=mode][value=pct]").checked=true; draw(); };
$("pctnum").oninput=()=>{ $("pct").value=$("pctnum").value; document.querySelector("input[name=mode][value=pct]").checked=true; draw(); };
$("thr").oninput=()=>{ $("thrnum").value=(+$("thr").value).toPrecision(4); document.querySelector("input[name=mode][value=thr]").checked=true; draw(); };
$("thrnum").oninput=()=>{ $("thr").value=$("thrnum").value; document.querySelector("input[name=mode][value=thr]").checked=true; draw(); };
document.querySelectorAll("input[name=mode]").forEach(e=>e.addEventListener("change",draw));
$("csv").onclick=()=>download(true); $("csvall").onclick=()=>download(false);
$("lochead").addEventListener("click",()=>{ $("loc").classList.toggle("big"); Plotly.Plots.resize("locplot"); });
document.body.className="v3d"; buildTree(); fillGraphSelect(); setThrRange(); draw();
</script></body></html>
"""
