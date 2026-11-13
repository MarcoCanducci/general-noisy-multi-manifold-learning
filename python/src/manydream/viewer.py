"""Standalone interactive viewer (one HTML file) for one or several fitted hierarchies.

    manydream.viewer.export_viewer("viewer.html", X, {"run A": H_a, "run B": H_b}, names=feature_names)

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
SCORE_HELP = {
    "logdens": "How likely the whole model (all structures together) makes a point at this position; "
               "low values mean the point lies where the model expects very few points.",
    "dist": "How far the point lies from its structure, measured across the structure (perpendicular "
            "to its local tangent directions), in the units of the data.",
    "d2": "The same perpendicular offset, measured in units of the noise the model fitted at that "
          "place, so that thick and thin parts of a structure can be compared.",
    "pval": "The probability, under the model's local noise, of a perpendicular offset at least this "
            "large; small values mean the offset is hard to explain as noise.",
    "rank": "The share of the neighbouring points along the structure (same node or within a few graph "
            "steps) that lie closer to the structure than this point; 1 means the farthest in its region.",
    "bg": "The model's probability that the point belongs to the uniform background and not to any structure.",
}


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


def _cloud_embedding(means: np.ndarray, edges: np.ndarray, P: np.ndarray, dim: int, k: int = 8,
                     max_landmarks: int = 1000, seed: int = 0):
    """Intrinsic coordinates of a structure measured through its data points.

    Nodes and points are joined in one neighbourhood graph (``k`` nearest neighbours plus the
    graph's own edges) and unfolded together by landmark Isomap. Unlike :func:`_graph_embedding`
    the layout does not depend on which edges the graph has, so points fill the structure and
    missing edges show as gaps in the graph rather than gaps in the picture.
    Returns the layout of the nodes and of the points.
    """
    from scipy.linalg import eigh
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import shortest_path
    from scipy.spatial import cKDTree

    K, m = len(means), int(min(max(dim, 2), 4))
    Z = np.vstack([means, P])
    n = len(Z)
    if n <= m + 1:
        return np.zeros((K, m)), np.zeros((len(P), m))
    kk = min(k, n - 1)
    dist, nb = cKDTree(Z).query(Z, k=kk + 1)
    i = np.repeat(np.arange(n), kk)
    j, w = nb[:, 1:].ravel(), dist[:, 1:].ravel()
    edges = np.asarray(edges, dtype=np.int64).reshape(-1, 2)
    i, j = np.r_[i, edges[:, 0]], np.r_[j, edges[:, 1]]
    w = np.r_[w, np.linalg.norm(means[edges[:, 0]] - means[edges[:, 1]], axis=1)] + 1e-12
    A = coo_matrix((w, (i, j)), shape=(n, n)).tocsr()
    A = A.maximum(A.T)
    land = np.arange(n) if n <= max_landmarks else np.sort(np.random.default_rng(seed).choice(n, max_landmarks, replace=False))
    Dl = shortest_path(A, method="D", directed=False, indices=land)
    finite = np.isfinite(Dl)
    Dl[~finite] = 1.5 * Dl[finite].max() if finite.any() else 1.0
    D2 = Dl[:, land] ** 2
    L = len(land)
    J = np.eye(L) - 1.0 / L
    vals, vecs = eigh(-0.5 * J @ D2 @ J, subset_by_index=[max(L - m, 0), L - 1])
    vals, vecs = np.maximum(vals[::-1], 1e-12), vecs[:, ::-1]
    emb = (-0.5 * (vecs / np.sqrt(vals)).T @ (Dl ** 2 - D2.mean(1, keepdims=True))).T
    out = np.zeros((n, m))
    out[:, :emb.shape[1]] = emb
    return out[:K], out[K:]


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


def _dimension_labels(H, n: int):
    """Estimated dimensionality of every point, before and after spatial smoothing.

    Read from the detection runs stored in ``H.detections_`` (a list of dicts with a fitted
    ``learner`` and the ``region`` of rows it was run on); where a finer level examined a point
    again its estimate replaces the coarser one. 0 marks points with no estimate: removed by
    the preprocessing (thinning, density filter) or with too small a neighbourhood.
    Returns ``None`` when the hierarchy holds no detection runs.
    """
    det = getattr(H, "detections_", None)
    if not det:
        return None
    hard, smooth = np.zeros(n, dtype=np.int8), np.zeros(n, dtype=np.int8)
    for d in sorted(det, key=lambda d: d.get("level", 0)):
        L = d["learner"]
        if getattr(L, "dim_result_", None) is None:
            continue
        rows = np.asarray(d["region"])[L.preprocess_.indices]
        ok = L.dim_result_.index > 0
        hard[rows[ok]] = L.dim_result_.index_hard[ok]
        smooth[rows[ok]] = L.dim_result_.index[ok]
    return hard, smooth


_STEP_LABELS = {"DensityThin": "thinned out (density subsampling)",
                "DensityFilter": "removed by the density filter (sparse neighbourhood)"}
USED = "used: on a structure's graph"


def _point_status(H, n: int):
    """What happened to every point during detection, and how many points each step set aside.

    Read from ``H.detections_``. Returns ``(codes, names, text)``: a status code per point (taken
    from the first-level runs; a later run on rows an earlier run removed replaces their status),
    the status names, and a table of counts for every detection run of every level. ``None``
    when no detection runs are stored. The fitted model still scores every point; the status
    says whether the point took part in building the structures.
    """
    det = getattr(H, "detections_", None)
    if not det:
        return None
    names = [USED]

    def code(name):
        if name not in names:
            names.append(name)
        return names.index(name)

    top = min(d.get("level", 0) for d in det)
    status = np.full(n, -1, dtype=np.int8)
    lines = []
    for d in sorted(det, key=lambda d: d.get("level", 0)):
        L, region = d["learner"], np.asarray(d["region"])
        if getattr(L, "dim_result_", None) is None:
            continue
        kept = np.asarray(L.preprocess_.indices)
        st = np.zeros(len(region), dtype=np.int8)
        rows = [("examined", len(region))]
        prev = np.arange(len(region))
        for name, idx in L.preprocess_.info.get("steps") or [("preprocessing", kept)]:
            gone = np.setdiff1d(prev, idx)
            if len(gone) or name in _STEP_LABELS:
                label = _STEP_LABELS.get(name, f"removed by {name}")
                st[gone] = code(label)
                rows.append((label, len(gone)))
            prev = np.asarray(idx)
        index = np.asarray(L.dim_result_.index)
        d_amb = L.X_.shape[1]
        dims = list(range(1, d_amb)) if L.dims is None else [int(j) for j in np.atleast_1d(L.dims)]
        und = index == 0
        skipped = ~und & ~np.isin(index, dims)
        on = np.zeros(len(kept), dtype=bool)
        for gs in getattr(L, "graphs_", {}).values():
            for g in gs:
                on[np.asarray(g.points, dtype=np.int64)] = True
        loose = ~und & ~skipped & ~on
        for mask, label in ((und, "dimensionality undetermined (too few neighbours)"),
                            (skipped, "estimated dimension not modelled"),
                            (loose, "not reached by any graph")):
            st[kept[mask]] = code(label)
            rows.append((label, int(mask.sum())))
        rows.append((USED, int((st == 0).sum())))
        if d.get("level", 0) == top:
            status[region] = st
        lines.append(f"level {d.get('level', '?')}, inside {d.get('parent', '?')}, r = {d.get('r', L.r):g}")
        lines += [f"  {c:>7,d}  {label}" for label, c in rows]
    if (status < 0).any():
        status[status < 0] = code("not examined")
    cnt = np.bincount(status, minlength=len(names))
    head = ["All points, after the first level:"] + [f"  {c:>7,d}  {nm}" for nm, c in zip(names, cnt) if c] + ["", "Every detection run:"]
    return status, names, "\n".join(head + lines)


def _rigid_fit(src: np.ndarray, dst: np.ndarray) -> np.ndarray:
    """``src`` turned (and mirrored if needed) and shifted onto ``dst`` in the least-squares sense."""
    a, b = src - src.mean(0), dst - dst.mean(0)
    U, _, Vt = np.linalg.svd(a.T @ b)
    return a @ (U @ Vt) + dst.mean(0)


def _topology(X, items, max_landmarks: int = 500) -> dict:
    """Layout of the whole mixture of graphs for the topology page.

    ``items`` holds, per graph that owns points: its structure id, node means, edges, the rows of
    its points and its own unfolding in up to three coordinates (nodes and points, through its
    data points; a curve or surface stays flat). All graphs and all points are also unfolded
    together in 3-D (one landmark Isomap through the whole point cloud), which gives every graph
    a position and orientation relative to the others.
    For both layouts the coordinates are stored relative to the graph's centre, so the page can
    move and turn graphs as rigid pieces. Pairs of graphs in contact (nodes of one within 1.5
    local edge lengths of a node of the other) are listed with the nodes involved.
    """
    from scipy.spatial import cKDTree

    n = len(X)
    out = {"islands": [], "links": [], "pj": "", "pi": "", "pf": "", "pisl": ""}
    if not items:
        return out
    off = np.cumsum([0] + [len(it["means"]) for it in items])
    P = np.vstack([it["means"] for it in items])
    E = np.vstack([np.asarray(it["edges"], dtype=np.int64).reshape(-1, 2) + off[k] for k, it in enumerate(items)])
    jn, jp = _cloud_embedding(P, E, X, 3, max_landmarks=max_landmarks)
    pj, pi = np.zeros((n, 3), dtype=np.float32), np.zeros((n, 3), dtype=np.float32)
    pisl = np.full(n, -1, dtype=np.int16)
    pf = np.zeros((n, 3), dtype=np.float32)
    el, placed = [], []
    for k, it in enumerate(items):
        rows = it["rows"]
        target = np.vstack([jn[off[k]:off[k + 1]], jp[rows]])
        c = target.mean(0)
        own = _rigid_fit(np.vstack([it["own_nodes"], it["own_points"]]), target)
        K = len(it["means"])
        placed.append(own)
        pj[rows], pi[rows], pisl[rows] = target[K:] - c, own[K:] - c, k
        out["islands"].append({"g": int(it["g"]), "c": [float(c[0]), float(c[1]), float(c[2])],
                               "nj": _b64(target[:K] - c, np.float32), "ni": _b64(own[:K] - c, np.float32)})
        e = np.asarray(it["edges"], dtype=np.int64).reshape(-1, 2)
        el.append(np.median(np.linalg.norm(it["means"][e[:, 0]] - it["means"][e[:, 1]], axis=1)) if len(e) else 0.0)
    el = np.array(el)
    el[el <= 0] = np.median(el[el > 0]) if (el > 0).any() else 1.0
    trees = [cKDTree(it["means"]) for it in items]
    seams = [[] for _ in items]                 # per graph: (own node ids, other graph, its node ids)
    for a in range(len(items)):
        for b in range(a + 1, len(items)):
            d, j = trees[b].query(items[a]["means"])
            ok = d < 1.5 * max(el[a], el[b])
            if ok.any():
                out["links"].append({"a": a, "b": b, "n": int(ok.sum()), "ia": _b64(np.flatnonzero(ok), np.int32),
                                     "ib": _b64(np.unique(j[ok]), np.int32)})
                d2, i2 = trees[a].query(items[b]["means"])
                ok2 = d2 < 1.5 * max(el[a], el[b])
                ia = np.r_[np.flatnonzero(ok), i2[ok2]]
                ib = np.r_[j[ok], np.flatnonzero(ok2)]
                seams[a].append((ia, b, ib))
                seams[b].append((ib, a, ia))
    # "flattened and joined": every graph keeps its own flat unfolding and is moved as a rigid
    # piece so that the nodes it shares a contact with meet those of its neighbours
    flat = _join_pieces([pl[:len(it["means"])] for pl, it in zip(placed, items)],
                        [jn[off[k]:off[k + 1]] for k in range(len(items))], seams)
    for k, it in enumerate(items):
        K = len(it["means"])
        moved = (placed[k] - flat[k]["cs"]) @ flat[k]["R"] + flat[k]["cd"]
        c = moved.mean(0)
        pf[it["rows"]] = moved[K:] - c
        out["islands"][k].update(cf=[float(c[0]), float(c[1]), float(c[2])], nf=_b64(moved[:K] - c, np.float32))
    out["seam_gap"] = {"islands": _seam_gap([pl[:len(it["means"])] for pl, it in zip(placed, items)], seams, el),
                       "joint": _seam_gap([jn[off[k]:off[k + 1]] for k in range(len(items))], seams, el),
                       "flat": _seam_gap([(placed[k][:len(it["means"])] - flat[k]["cs"]) @ flat[k]["R"] + flat[k]["cd"]
                                          for k, it in enumerate(items)], seams, el)}
    out.update(pj=_b64(pj, np.float32), pi=_b64(pi, np.float32), pf=_b64(pf, np.float32), pisl=_b64(pisl, np.int16))
    return out


def _seam_gap(nodes, seams, el) -> float:
    """Mean distance between nodes in contact, in units of the local edge length."""
    g = [np.linalg.norm(nodes[a][ia] - nodes[b][ib], axis=1) / max(el[a], el[b])
         for a, lst in enumerate(seams) for ia, b, ib in lst if a < b]
    return float(np.concatenate(g).mean()) if g else 0.0


def _join_pieces(nodes, anchors, seams, anchor_weight: float = 0.05, sweeps: int = 80) -> list:
    """Rigid motions that make flat pieces meet along their contacts.

    ``nodes[k]`` are the node positions of piece ``k`` (already placed roughly), ``anchors[k]``
    where the joint unfolding puts them, ``seams[k]`` its contacts ``(own ids, other piece, its
    ids)``. Each piece in turn is moved (rotation and shift, no deformation) to bring its contact
    nodes onto those of its neighbours, with a weak pull towards its anchors that keeps the
    overall arrangement and settles pieces without contacts. Returns per piece ``R``, ``cs``,
    ``cd`` with ``moved = (x - cs) @ R + cd``.
    """
    cur = [a.copy() for a in nodes]
    T = [{"R": np.eye(3), "cs": np.zeros(3), "cd": np.zeros(3)} for _ in nodes]
    for _ in range(sweeps):
        shift = 0.0
        for k, base in enumerate(nodes):
            src = [base] + [base[ia] for ia, _, _ in seams[k]]
            dst = [anchors[k]] + [cur[b][ib] for _, b, ib in seams[k]]
            w = np.concatenate([np.full(len(base), anchor_weight)] + [np.ones(len(ia)) for ia, _, _ in seams[k]])
            A, B = np.vstack(src), np.vstack(dst)
            cs, cd = (w[:, None] * A).sum(0) / w.sum(), (w[:, None] * B).sum(0) / w.sum()
            U, _, Vt = np.linalg.svd((A - cs).T @ (w[:, None] * (B - cd)))
            R = U @ Vt
            new = (base - cs) @ R + cd
            shift = max(shift, float(np.abs(new - cur[k]).max()))
            cur[k], T[k] = new, {"R": R, "cs": cs, "cd": cd}
        if shift < 1e-5:
            break
    return T


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
    ppos2 = np.zeros((len(X), 4), dtype=np.float32)
    graph_of = np.array([s["graph_of"] for s in structures] + [-1])
    gid = graph_of[sid]
    models = {node_id[k]: None for k in node_id if not isinstance(k, tuple)}

    def collect(node):
        for c in node.children:
            models[node_id[id(c)]] = c.model
            collect(c)

    collect(H.root_)
    items = []
    for g_, model in models.items():
        rows = np.flatnonzero(gid == g_)
        if len(rows) == 0:
            continue
        emb = _graph_embedding(model.means_, model.graph.edges, int(model.graph.dim))[0]
        ppos[rows, :emb.shape[1]] = _point_layout(X[rows], model.means_, emb, model.graph.edges,
                                                  np.asarray(res["node"])[rows], int(model.graph.dim))
        if int(model.graph.dim) >= 2:
            e2, p2 = _cloud_embedding(model.means_, model.graph.edges, X[rows], int(model.graph.dim))
            structures[g_]["emb2"] = _b64(e2, np.float32)
            ppos2[rows, :p2.shape[1]] = p2
        else:
            e2, p2 = _cloud_embedding(model.means_, model.graph.edges, X[rows], 2)
        pad = lambda a: np.c_[a, np.zeros((len(a), 3 - a.shape[1]))] if a.shape[1] < 3 else a[:, :3]  # noqa: E731
        items.append({"g": g_, "means": model.means_, "edges": model.graph.edges, "rows": rows,
                      "own_nodes": pad(e2), "own_points": pad(p2)})
    counts = np.bincount(sid[sid >= 0], minlength=len(structures))
    for s in structures:
        s["points"] = int(counts[s["id"]])
    idx = _dimension_labels(H, len(X))
    stat = _point_status(H, len(X))
    return {
        "status": _b64(stat[0], np.int8) if stat else "",
        "status_names": stat[1] if stat else [],
        "discarded": stat[2] if stat else "No detection runs are stored for this version.",
        "idxh": _b64(idx[0], np.int8) if idx else "",
        "idxs": _b64(idx[1], np.int8) if idx else "",
        "structures": structures,
        "sid": _b64(sid, np.int16),
        "node": _b64(res["node"], np.int32),
        "toward": _b64(res["toward"], np.int32),
        "frac": _b64(res["fraction"], np.float32),
        "ppos": _b64(ppos, np.float32),
        "ppos2": _b64(ppos2, np.float32),
        "topo": _topology(X, items),
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
            "scores": [{"key": k, "label": l, "low": low, "help": SCORE_HELP.get(k, "")} for k, l, low in SCORES],
            "versions": [{"name": str(name), **_version_payload(H, X, hops, min_region)} for name, H in versions.items()]}


def export_viewer(path, X, versions, names=None, title: str = "MANY-DREAM viewer", plotlyjs: str = "inline",
                  hops: int = 2, min_region: int = 30, extra=None) -> Path:
    """Write the standalone viewer page.

    Parameters
    ----------
    path : file to write (``.html``)
    X : (n, d) data
    versions : dict name -> fitted :class:`manydream.HierarchicalAGTM` (or a single hierarchy)
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
.info{position:relative;display:inline-block;width:15px;height:15px;line-height:15px;text-align:center;border-radius:50%;border:1px solid var(--muted);color:var(--muted);font-size:10px;font-weight:600;cursor:help;margin-left:4px}
.infobox{display:none;position:fixed;left:16px;top:44px;width:310px;max-height:calc(100vh - 80px);overflow:auto;z-index:20;background:#ffffff;color:#0b0b0b;border:1px solid var(--line);border-radius:6px;box-shadow:0 4px 16px rgba(0,0,0,.18);padding:10px 12px;text-align:left;font-weight:400;font-size:11.5px;line-height:1.45;cursor:default}
.info:hover .infobox,.info:focus .infobox{display:block}.infobox b{display:block;margin-top:6px}.infobox b:first-child{margin-top:0}
#catlegend{display:flex;flex-wrap:wrap;gap:4px 10px;margin-top:4px;font-size:11px}#catlegend span{display:flex;align-items:center;gap:4px}
.viewonly{display:none}body.v3d .only3d{display:block}body.vgraph .onlygraph{display:block}body.vtopo .onlytopo{display:block}
#topo{display:none;flex:1 1 0;min-height:0;width:100%;background:#ffffff;cursor:grab;touch-action:none}body.vtopo #topo{display:block}body.vtopo #plot{display:none}
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
  <label for="colour">Colour points by</label><select id="colour"></select><div id="catlegend"></div>
  <div class="row"><div><label for="psize">Point size</label><input id="psize" type="range" min="1" max="6" step="0.5" value="2"></div>
  <div><label for="popac">Opacity</label><input id="popac" type="range" min="0.1" max="1" step="0.05" value="0.7"></div></div>
  <label class="check"><input type="checkbox" id="showpts" checked> points</label>
  <label class="check"><input type="checkbox" id="showgraphs" checked> graphs of the selected structures</label>
  <div class="viewonly onlytopo">
  <h2>Topology</h2>
  <label for="tmode">Layout</label><select id="tmode"><option value="isl">islands: each structure in its own unfolding</option><option value="flat">flattened and joined: flat pieces meeting at their contacts</option><option value="joint">joint: everything unfolded together</option></select>
  <label for="tsep">Pull the structures apart</label><input id="tsep" type="range" min="0" max="1" step="0.01" value="1">
  <label for="tlink">Hide links with fewer contacts than <span id="tlinkv">1</span></label><input id="tlink" type="range" min="1" max="50" step="1" value="1">
  <label for="tnode">Colour nodes by</label><select id="tnode"></select>
  <label class="check"><input type="checkbox" id="tpts" checked> data points</label>
  <label class="check"><input type="checkbox" id="tlinks" checked> links between structures in contact</label>
  <label class="check"><input type="checkbox" id="tlabels" checked> labels and outlines</label>
  <label class="check"><input type="checkbox" id="tpersp" checked> perspective (near things larger)</label>
  <div class="row tight" style="margin-top:6px"><button id="treset">Reset arrangement</button><button id="tfit">Fit to window</button><button id="tpng">Save picture</button></div>
  <p class="hint">The picture is three-dimensional. Drag the background to rotate the whole frame: the structures keep the arrangement you gave them and turn together. Drag a structure to move it (in the plane of the screen), Shift-drag to turn it on itself, press F over it to mirror it. Right-drag (or Ctrl-drag) to pan, use the wheel to zoom, double-click a structure to open its graph page. Both layouts are unfolded through the original data points: "islands" draws every structure in its own unfolding as on the graph page (curves and surfaces stay flat sheets), placed and oriented as in the joint unfolding; "flattened and joined" keeps every structure flat in the same way but moves each one, as a rigid piece, so that the nodes it has in contact with its neighbours meet (nothing is bent or stretched); "joint" unfolds all structures and all points together in 3-D, so structures that run through each other overlap and bend. A link joins two structures that have nodes within 1.5 edge lengths of each other; its width is the number of such nodes. Your arrangement is remembered in this browser.</p>
  </div>
  <div class="viewonly onlygraph">
  <h2>Graph view</h2>
  <label for="gstruct">Structure</label><select id="gstruct"></select>
  <label for="gcolour">Colour nodes by</label><select id="gcolour"></select>
  <label for="gedge">Colour edges by</label><select id="gedge"><option value="none">nothing (grey)</option><option value="curv">curvature along the edge</option><option value="len">edge length in data space</option><option value="stretch">distortion of the layout</option></select>
  <div id="gprojwrap" style="display:none"><label for="gproj">Projection (three of the four intrinsic coordinates)</label><select id="gproj"><option value="123">coordinates 1, 2, 3</option><option value="124">coordinates 1, 2, 4</option><option value="134">coordinates 1, 3, 4</option><option value="234">coordinates 2, 3, 4</option></select></div>
  <label class="check"><input type="checkbox" id="gsize" checked> node size by number of points</label>
  <label for="gy">Vertical axis (curves)</label><select id="gy"></select>
  <label for="glayout">Layout (surfaces and higher)</label><select id="glayout"><option value="cloud">unfolded through the data points</option><option value="graph">unfolded along the graph edges</option></select>
  <label class="check"><input type="checkbox" id="gpts" checked> points of the structure</label>
  <label class="check"><input type="checkbox" id="gfill" checked> fill the structure (points at their own position)</label>
  <label class="check"><input type="checkbox" id="gloc" checked> locator: where the structure sits in the data</label>
  <p class="hint">Each structure is drawn in its own dimension, unfolded along its graph (geodesic distances along the edges): a curve as a line with a profile above it, a surface in the plane, a 3-D structure in space, a 4-D one through 3-D projections. Two layouts: "through the data points" unfolds nodes and original data points together (distances through the point cloud), so points fill the structure and missing edges show as gaps in the graph; "along the graph edges" unfolds the graph alone and places each point by its offset from the nearest node, which pulls apart branches that the graph does not connect.</p>
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
  <div id="xinfo" class="hint"></div>
  <h2>Outliers</h2>
  <label for="score">Score <span class="info" tabindex="0">i<span class="infobox" id="scorehelp"></span></span></label><select id="score"></select>
  <label class="check"><input type="radio" name="mode" value="pct" checked> by percentile</label>
  <div class="row"><input id="pct" type="range" min="0.05" max="20" step="0.05" value="1"><input id="pctnum" type="number" min="0.01" max="50" step="0.05" value="1" style="max-width:90px"></div>
  <label class="check"><input type="radio" name="mode" value="thr"> by threshold</label>
  <div class="row"><input id="thr" type="range" step="any"><input id="thrnum" type="number" step="any" style="max-width:110px"></div>
  <label for="scope">Applied</label><select id="scope"><option value="shown">over all shown points</option><option value="each">within each shown structure</option></select>
  <label class="check"><input type="checkbox" id="hl" checked> highlight outliers</label>
  <label class="check"><input type="checkbox" id="only"> show outliers only</label>
  <div id="outinfo"></div>
  <div class="row tight" style="margin-top:6px"><button id="csv">Download outliers (CSV)</button><button id="csvall">Download shown points (CSV)</button></div>
  <details><summary>Points set aside at each step</summary><pre id="discarded"></pre></details>
  <details><summary>Model summary of this version</summary><pre id="summary"></pre></details>
  <h2>Table of all criteria</h2>
  <p class="hint">One line for every point of the data set (not only the shown or selected ones), with its status in the analysis and the criteria ticked here.</p>
  <div id="tabcrit"></div>
  <label class="check"><input type="checkbox" id="tabstruct" checked> structure, node and estimated dimensionality</label>
  <label class="check"><input type="checkbox" id="tabsel" checked> current selection (outlier yes/no)</label>
  <label class="check"><input type="checkbox" id="tabx"> coordinates</label>
  <label class="check"><input type="checkbox" id="tabextra"> extra labels</label>
  <p class="hint" style="margin-top:6px">Write NaN in the criteria of points with status:</p>
  <div id="tabnan"></div>
  <div class="row tight" style="margin-top:6px"><button id="tabdl">Download table, all points (CSV)</button></div>
</div>
<div id="main"><div id="tabs"><button class="tab active" data-view="3d">3-D view</button><button class="tab" data-view="graph">Graph view (one structure, unfolded)</button><button class="tab" data-view="topo">Topology (all structures)</button></div>
<div id="plot"></div><canvas id="topo"></canvas><div id="status"></div>
<div id="loc" class="viewonly onlygraph"><div id="lochead" title="click to enlarge or shrink"><span id="loctitle"></span><span>&#x2922;</span></div><div id="locplot"></div></div></div>
<script>
const DATA = /*__DATA__*/;
const PAL = ["#2a78d6","#eb6834","#1baf7a","#eda100","#e87ba4","#008300","#6250d6","#e34948"];
const DIMCOL = ["#b9b8b2","#e34948","#2a78d6","#1baf7a","#eda100","#6250d6","#e87ba4","#008300","#52514e"];
function dec(b64, T){ if(!b64) return new T(0); const s=atob(b64), u=new Uint8Array(s.length); for(let i=0;i<s.length;i++)u[i]=s.charCodeAt(i); return new T(u.buffer); }
function colourOf(i){ if(i<PAL.length) return PAL[i]; const h=(i*137.508)%360; return `hsl(${h.toFixed(0)},62%,${i%2?44:56}%)`; }
const N=DATA.n, D=DATA.d, X=dec(DATA.X, Float32Array);
const V = DATA.versions.map(v => ({name:v.name, summary:v.summary, sid:dec(v.sid,Int16Array), node:dec(v.node,Int32Array),
  toward:dec(v.toward,Int32Array), frac:dec(v.frac,Float32Array), ppos: v.ppos ? dec(v.ppos,Float32Array) : null, ppos2: v.ppos2 ? dec(v.ppos2,Float32Array) : null,
  status: v.status ? dec(v.status,Int8Array) : null, status_names: v.status_names||[], discarded: v.discarded||"",
  topo: v.topo && v.topo.islands && v.topo.islands.length ? {pj:dec(v.topo.pj,Float32Array), pi:dec(v.topo.pi,Float32Array), pf: v.topo.pf ? dec(v.topo.pf,Float32Array) : null, gap:v.topo.seam_gap||null, pisl:dec(v.topo.pisl,Int16Array),
    isl:v.topo.islands.map(t=>({g:t.g, c:t.c, nj:dec(t.nj,Float32Array), ni:dec(t.ni,Float32Array), nf: t.nf ? dec(t.nf,Float32Array) : null, cf:t.cf||t.c, d:[0,0,0], R:[1,0,0,0,1,0,0,0,1]})),
    links:v.topo.links.map(l=>({a:l.a,b:l.b,n:l.n,ia:dec(l.ia,Int32Array),ib:dec(l.ib,Int32Array)}))} : null,
  idxh: v.idxh ? dec(v.idxh,Int8Array) : null, idxs: v.idxs ? dec(v.idxs,Int8Array) : null,
  scores:Object.fromEntries(Object.entries(v.scores).map(([k,b])=>[k,dec(b,Float32Array)])),
  S:v.structures.map(s=>({...s, nodes:dec(s.nodes,Float32Array), edges:dec(s.edges,Int32Array), emb:dec(s.emb,Float32Array), emb2: s.emb2 ? dec(s.emb2,Float32Array) : null, arc:dec(s.arc,Float32Array), curv:dec(s.curv,Float32Array), zeta:dec(s.zeta,Float32Array), degen:dec(s.degen,Float32Array), noise:dec(s.noise,Float32Array), cov:dec(s.cov,Float32Array), on:s.kind!=="group", colour:colourOf(s.id)}))}));
let view = "3d";
const TRI = D*(D+1)/2;
let lastThree = null;
const EXTRA = (DATA.extra||[]).map(e=>({name:e.name, kind:e.kind, cats:e.cats||[], values: e.kind==="num" ? dec(e.values,Float32Array) : dec(e.values,Int32Array)}));
function extraOf(key){ return key && key.startsWith("x:") ? EXTRA[+key.slice(2)] : null; }
function refreshExtraOptions(){
  const keep=(id,items)=>{ const el=$(id), old=el.value; fill(el, items); if(items.some(([v])=>v===old)) el.value=old; };
  const xs=EXTRA.map((e,i)=>["x:"+i, e.name+(e.kind==="cat"?" (categories)":"")]), xn=xs.filter((_,i)=>EXTRA[i].kind==="num");
  keep("colour", [["structure","structure"],["dim","dimension of the structure"],["idxh","estimated dimensionality, before smoothing"],["idxs","estimated dimensionality, after smoothing"]].concat(DATA.scores.map(s=>[s.key,s.label]), xs));
  keep("gcolour", [["count","points per node"],["median","median of the outlier score"],["nout","outliers per node"],["deg","node degree"],["curv","curvature at the node"],["noise","noise scale of the node"],["zeta","covariance scale zeta"],["degen","degenerate components"]].concat(DATA.names.map((n,i)=>["pc:"+i,"mean "+n]), xn.map(([v,t])=>[v,"mean "+t])));
  keep("tnode", [["structure","the structure's colour"]].concat([...$("gcolour").options].map(o=>[o.value,o.text])));
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
  $("summary").textContent=v.summary; $("discarded").textContent=v.discarded;
  const nb=$("tabnan"); nb.innerHTML=""; v.status_names.forEach((t,k)=>{ if(k===0) return; const l=document.createElement("label"); l.className="check"; const c=document.createElement("input"); c.type="checkbox"; c.dataset.code=k; c.checked=!t.startsWith("thinned"); l.append(c," "+t); nb.appendChild(l); });
  if(!v.status) nb.textContent="(no status stored for this version)";
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
function draw(){ buildDimToggles(); catLegend(); if(view==="graph") drawGraph(); else if(view==="topo") drawTopo(true); else draw3d(); }
function catLegend(){ const cm=$("colour").value, e=extraOf(cm), box=$("catlegend"); box.innerHTML="";
  if(cm==="idxh"||cm==="idxs"||cm==="dim"){ const a=V[cur][cm], cnt=new Array(9).fill(0);
    if(cm==="dim"){ for(let i=0;i<N;i++) cnt[Math.min(V[cur].S[V[cur].sid[i]]?.dim||0,8)]++; } else if(a){ for(let i=0;i<N;i++) cnt[Math.min(a[i],8)]++; } else { box.textContent="No dimensionality estimates are stored for this version."; return; }
    cnt.forEach((c,d)=>{ if(!c) return; const sp=document.createElement("span"), sw=document.createElement("span"); sw.className="sw"; sw.style.background=DIMCOL[d];
      sp.append(sw, document.createTextNode((d? d+"-D" : (cm==="dim"?"none":"not estimated"))+" ("+c.toLocaleString()+")")); box.appendChild(sp); }); return; }
  if(!e||e.kind!=="cat") return;
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
function graphPointColour(v,idx,nr){
  // colours of the structure's points on the graph page, following "Colour points by"; null = plain
  const cm=$("colour").value, e=extraOf(cm);
  if(cm==="idxh"||cm==="idxs") return v[cm] ? {f:i=>DIMCOL[Math.min(v[cm][i],8)], marker:{}} : null;
  if(e && e.kind==="cat") return {f:i=> e.values[i]<0 ? "#b9b8b2" : colourOf(e.values[i]), marker:{}};
  const vals = e ? e.values : v.scores[cm]; if(!vals) return null;
  // same scale as the nodes when they show the mean of the same property
  const rg = $("gcolour").value===cm ? nr : ramp(idx.map(i=>vals[i]));
  return {f:i=>vals[i], marker:{colorscale:RAMPSCALE,cmin:rg.lo,cmax:rg.hi}};
}
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
    const px=[],py=[],cd=[],ox=[],oy=[],ocd=[],pcol=[],gpc=graphPointColour(v,idx,nr);
    const yk=$("gy").value, yf = yk==="score" ? (i=>sc[i]) : pointValue(yk), ylab = yk==="score" ? meta.label : (yk.startsWith("pc:") ? DATA.names[+yk.slice(3)] : extraOf(yk).name);
    idx.forEach(i=>{ const a=v.node[i], b=v.toward[i]; const x=g.arc[a]+Math.sign(g.arc[b]-g.arc[a])*v.frac[i]*elen(a,b), y=yf(i);
      if(!Number.isFinite(y)) return;
      if(r.out[i]){ ox.push(x);oy.push(y);ocd.push([i,sc[i]]); if(hl||only) return; } if(only) return; px.push(x);py.push(y);cd.push([i,sc[i]]); if(gpc) pcol.push(gpc.f(i)); });
    if(showPts) traces.push({type:"scattergl",mode:"markers",x:px,y:py,customdata:cd,hovertemplate:ht,marker:{...(gpc?gpc.marker:{}),size:5,color:gpc?pcol:s.colour,opacity:0.7},showlegend:false});
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
    const cloud = $("glayout").value==="cloud" && g.emb2 && v.ppos2, EMB = cloud ? g.emb2 : g.emb, PP = cloud ? v.ppos2 : v.ppos;
    const nd=three?3:2, pos=(k,a)=>EMB[k*m+ax[a]], T=three?"scatter3d":"scattergl";
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
    const jit=0.25*(median(med)||0.01), fill=cloud || ($("gfill").checked && v.ppos);
    const px=[],py=[],pz=[],cd=[],ox=[],oy=[],oz=[],ocd=[],pcol=[], pidx = graphPointColour(v,idx,nr);
    idx.forEach(i=>{ const a=v.node[i], b=v.toward[i], f=v.frac[i], c=[];
      for(let q=0;q<nd;q++) c.push(fill ? PP[i*4+ax[q]] : pos(a,q)+f*(pos(b,q)-pos(a,q))+jit*(hash01(i,q+1)-0.5)*1.6);
      if(r.out[i]){ ox.push(c[0]);oy.push(c[1]);if(three)oz.push(c[2]);ocd.push([i,sc[i]]); if(hl||only) return; } if(only) return; px.push(c[0]);py.push(c[1]);if(three)pz.push(c[2]);cd.push([i,sc[i]]); if(pidx) pcol.push(pidx.f(i)); });
    if(showPts && px.length){ const t={type:T,mode:"markers",x:px,y:py,customdata:cd,hovertemplate:ht,marker:{...(pidx?pidx.marker:{}),size:three?(pidx?2.2:1.6):(pidx?4:3),color:pidx?pcol:"#b9b8b2",opacity:pidx?0.85:(three?0.5:0.6),line:{width:0}},showlegend:false}; if(three) t.z=pz; traces.push(t); }
    const nx=[],ny=[],nz=[],ncd=[],sz=[]; for(let k=0;k<K;k++){ nx.push(pos(k,0)); ny.push(pos(k,1)); if(three) nz.push(pos(k,2)); ncd.push([k,np.cnt[k],np.nout[k],np.val[k]]); sz.push(three? nsize(k)*0.55 : nsize(k)); }
    const nt={type:T,mode:"markers",x:nx,y:ny,customdata:ncd,hovertemplate:nht,marker:nodeMarker(sz),showlegend:false}; if(three){ nt.z=nz; nt.marker.line={width:0}; } traces.push(nt);
    if(ox.length){ const t={type:T,mode:"markers",x:ox,y:oy,customdata:ocd,hovertemplate:ht,marker:{size:three?4.5:9,color:"#e34948",symbol:"diamond",line:{color:"#fff",width:0.5}},showlegend:false}; if(three) t.z=oz; traces.push(t); }
    if(three){ const axn=a=>({title:{text:"intrinsic coordinate "+(ax[a]+1)},backgroundcolor:"#ffffff",gridcolor:"#eeede8",zerolinecolor:"#eeede8",color:"#52514e",showticklabels:false});
      layout={scene:{aspectmode:"data",xaxis:axn(0),yaxis:axn(1),zaxis:axn(2)},margin:{l:0,r:0,t:0,b:0}}; }
    else layout={xaxis:{visible:false},yaxis:{visible:false,scaleanchor:"x"},margin:{l:10,r:90,t:10,b:10},dragmode:"pan"};
  }
  Plotly.react("plot", traces, Object.assign({paper_bgcolor:"#fcfcfb",plot_bgcolor:"#ffffff",uirevision:"g"+cur+"-"+sidSel+"-"+$("gproj").value+$("glayout").value},layout), {displaylogo:false,responsive:true,scrollZoom:true});
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

// ---------------------------------------------------------------- topology page (3-D, canvas)
const I3=()=>[1,0,0,0,1,0,0,0,1];
const mm=(A,B)=>{ const C=new Array(9); for(let i=0;i<3;i++) for(let j=0;j<3;j++) C[3*i+j]=A[3*i]*B[j]+A[3*i+1]*B[3+j]+A[3*i+2]*B[6+j]; return C; };
const mt=A=>[A[0],A[3],A[6],A[1],A[4],A[7],A[2],A[5],A[8]];
const mv=(A,x)=>[A[0]*x[0]+A[1]*x[1]+A[2]*x[2], A[3]*x[0]+A[4]*x[1]+A[5]*x[2], A[6]*x[0]+A[7]*x[1]+A[8]*x[2]];
function rotXY(ax,ay){ // rotation by ax about the screen x axis and ay about the screen y axis
  const cx=Math.cos(ax), sx=Math.sin(ax), cy=Math.cos(ay), sy=Math.sin(ay); return mm([1,0,0,0,cx,-sx,0,sx,cx],[cy,0,sy,0,1,0,-sy,0,cy]); }
function orthonormal(A){ // keep accumulated rotations rigid
  let x=[A[0],A[3],A[6]], y=[A[1],A[4],A[7]]; const nrm=u=>{ const l=Math.hypot(u[0],u[1],u[2])||1; return [u[0]/l,u[1]/l,u[2]/l]; }; x=nrm(x);
  const d=x[0]*y[0]+x[1]*y[1]+x[2]*y[2]; y=nrm([y[0]-d*x[0],y[1]-d*x[1],y[2]-d*x[2]]); let z=[x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1]-x[1]*y[0]];
  const det=A[0]*(A[4]*A[8]-A[5]*A[7])-A[1]*(A[3]*A[8]-A[5]*A[6])+A[2]*(A[3]*A[7]-A[4]*A[6]); if(det<0) z=[-z[0],-z[1],-z[2]];
  return [x[0],y[0],z[0],x[1],y[1],z[1],x[2],y[2],z[2]]; }
const TV={scale:1,tx:0,ty:0,R:I3(),fitted:""}; let topoHover=-1, topoDrag=null, topoCache=null, topoRAF=0;
function topoKey(){ return "manydream-topo3|"+document.title+"|"+V[cur].name+"|"+$("tmode").value; }
function topoSave(){ try{ localStorage.setItem(topoKey(), JSON.stringify({isl:V[cur].topo.isl.map(t=>[t.d,t.R]), view:TV.R})); }catch(e){} }
function topoLoad(){ const T=V[cur].topo; if(!T) return; let a=null; try{ a=JSON.parse(localStorage.getItem(topoKey())||"null"); }catch(e){}
  const ok=a && a.isl && a.isl.length===T.isl.length; T.isl.forEach((t,k)=>{ t.d = ok ? a.isl[k][0].slice() : [0,0,0]; t.R = ok ? a.isl[k][1].slice() : I3(); }); TV.R = ok && a.view ? a.view.slice() : rotXY(-0.35,0.45); }
function topoGeometry(){
  // radii, separated centres and link anchors of the current version and layout (cached)
  const v=V[cur], T=v.topo, mode=$("tmode").value, key=cur+"|"+mode; if(T.geo && T.geo.key===key) return T.geo;
  const K=T.isl.length, R=new Array(K), flat = mode==="flat" && T.pf && T.isl.every(t=>t.nf), loc=T.isl.map(t=> mode==="joint" ? t.nj : flat ? t.nf : t.ni), pl = mode==="joint" ? T.pj : flat ? T.pf : T.pi, cen=T.isl.map(t=> flat ? t.cf : t.c);
  const far=T.isl.map(()=>[]); for(let k=0;k<K;k++){ const L=loc[k]; for(let i=0;i<L.length/3;i++) far[k].push(Math.hypot(L[3*i],L[3*i+1],L[3*i+2])); }
  const step=Math.max(1,Math.floor(N/6000)); for(let i=0;i<N;i+=step){ const k=T.pisl[i]; if(k>=0) far[k].push(Math.hypot(pl[3*i],pl[3*i+1],pl[3*i+2])); }
  for(let k=0;k<K;k++) R[k]=quantile(far[k],0.97)||0.1;
  const medR=median(R), Rm=R.map(x=>Math.max(x,0.3*medR)), gap=0.12*medR, C=cen.map(c=>c.slice()), C0=C.map(c=>c.slice());
  for(let it=0;it<600;it++){ let moved=false;
    for(let a=0;a<K;a++) for(let b=a+1;b<K;b++){ let d=[C[b][0]-C[a][0],C[b][1]-C[a][1],C[b][2]-C[a][2]], L=Math.hypot(d[0],d[1],d[2]); if(L<1e-9){ d=[1e-3*(b+1),1e-3,1e-3*a]; L=Math.hypot(d[0],d[1],d[2]); }
      const need=Rm[a]+Rm[b]+gap; if(L<need){ const sh=(need-L)/L, wa=Rm[b]/(Rm[a]+Rm[b]); for(let q=0;q<3;q++){ C[a][q]-=d[q]*sh*wa; C[b][q]+=d[q]*sh*(1-wa); } moved=true; } }
    if(!moved) break; if(it<400) for(let k=0;k<K;k++) for(let q=0;q<3;q++) C[k][q]+=0.01*(C0[k][q]-C[k][q]); }
  const anchor=(L,ids)=>{ const m=[0,0,0]; ids.forEach(i=>{ m[0]+=L[3*i]; m[1]+=L[3*i+1]; m[2]+=L[3*i+2]; }); return m.map(x=>x/ids.length); };
  const links=T.links.map(l=>({a:l.a,b:l.b,n:l.n,pa:anchor(loc[l.a],l.ia),pb:anchor(loc[l.b],l.ib)}));
  const rows=T.isl.map(()=>[]); for(let i=0;i<N;i++){ const k=T.pisl[i]; if(k>=0) rows[k].push(i); }
  T.geo={key,R:Rm,sep:C,cen,loc,pl,links,rows,maxn:Math.max(1,...links.map(l=>l.n))}; return T.geo;
}
function topoCentre(k,G){ const t=V[cur].topo.isl[k], c=G.cen[k], f=+$("tsep").value; return [0,1,2].map(q=>c[q]+f*(G.sep[k][q]-c[q])+t.d[q]); }
function topoPersp(G){ // perspective about the middle of the window: size factor for a depth z in view coordinates
  const cv=$("topo"); if(!$("tpersp").checked) return {on:false,f:()=>1,cx:0,cy:0};
  let z0=0,n=0,ext=0; const cs=[]; V[cur].topo.isl.forEach((t,k)=>{ if(!topoVisible(k)) return; const c=mv(TV.R,topoCentre(k,G)); cs.push([c[2],G.R[k]]); z0+=c[2]; n++; });
  z0 = n ? z0/n : 0; cs.forEach(c=>{ ext=Math.max(ext,Math.abs(c[0]-z0)+c[1]); });
  const dist=2.2*Math.max(ext,1e-9); return {on:true, cx:cv.clientWidth/2, cy:cv.clientHeight/2, f:z=>dist/Math.max(dist-(z-z0),0.15*dist)}; }
function topoScreen(c,P){ const f=P.f(c[2]), x=TV.tx+TV.scale*c[0], y=TV.ty-TV.scale*c[1]; return P.on ? [P.cx+(x-P.cx)*f, P.cy+(y-P.cy)*f, f] : [x,y,1]; }
function topoVisible(k){ const v=V[cur], g=v.topo.isl[k].g; return v.S.some(s=>s.on && s.kind!=="group" && s.graph_of===g); }
function topoFit(){ const cv=$("topo"), v=V[cur]; if(!v.topo) return; const G=topoGeometry();
  let x0=Infinity,x1=-Infinity,y0=Infinity,y1=-Infinity; v.topo.isl.forEach((t,k)=>{ if(!topoVisible(k)) return; const p=mv(TV.R,topoCentre(k,G)), r=G.R[k]*1.15;
    x0=Math.min(x0,p[0]-r); x1=Math.max(x1,p[0]+r); y0=Math.min(y0,p[1]-r); y1=Math.max(y1,p[1]+r); });
  if(!Number.isFinite(x0)) return; const w=cv.clientWidth, h=cv.clientHeight; TV.scale=0.94*Math.min(w/(x1-x0), h/(y1-y0)); TV.tx=w/2-TV.scale*(x0+x1)/2; TV.ty=h/2+TV.scale*(y0+y1)/2; }
function topoHit(mx,my){ const v=V[cur]; if(!v.topo) return -1; const G=topoGeometry(), P=topoPersp(G); let best=-1, bz=-Infinity;
  v.topo.isl.forEach((t,k)=>{ if(!topoVisible(k)) return; const p=mv(TV.R,topoCentre(k,G)), q=topoScreen(p,P), sx=q[0], sy=q[1];
    if(Math.hypot(mx-sx,my-sy)<=G.R[k]*1.05*TV.scale*q[2] && (best<0 || G.R[k]<G.R[best]*0.999 || (Math.abs(G.R[k]-G.R[best])<1e-9 && p[2]>bz))){ best=k; bz=p[2]; } }); return best; }
function topoColours(r){
  // colour index of every point and the palette, following "Colour points by"
  const v=V[cur], cm=$("colour").value, e=extraOf(cm), ci=new Uint16Array(N), pal=[], idOf=new Map(); let bar=null;
  const id=c=>{ let q=idOf.get(c); if(q===undefined){ q=pal.length; pal.push(c); idOf.set(c,q); } return q; };
  if(cm==="structure") for(let i=0;i<N;i++) ci[i]=id(v.S[v.sid[i]] ? v.S[v.sid[i]].colour : "#b9b8b2");
  else if(cm==="dim") for(let i=0;i<N;i++) ci[i]=id(DIMCOL[Math.min(v.S[v.sid[i]] ? v.S[v.sid[i]].dim : 0,8)]);
  else if(cm==="idxh"||cm==="idxs") for(let i=0;i<N;i++) ci[i]=id(v[cm] ? DIMCOL[Math.min(v[cm][i],8)] : "#b9b8b2");
  else if(e && e.kind==="cat") for(let i=0;i<N;i++) ci[i]=id(e.values[i]<0 ? "#b9b8b2" : colourOf(e.values[i]));
  else { const vals = e ? e.values : v.scores[cm], f=[]; for(let i=0;i<N;i++) if(r.shown[i] && Number.isFinite(vals[i])) f.push(vals[i]); const rg=ramp(f);
    for(let i=0;i<N;i++){ const b=rg.bin(vals[i]); ci[i]=id(b<0 ? "#d6d5cf" : RAMP[b]); } bar={title:"points: "+(e ? e.name : DATA.scores.find(q=>q.key===cm).label), lo:rg.lo, hi:rg.hi}; }
  return {ci,pal,bar};
}
function drawTopo(recompute){
  const cv=$("topo"), v=V[cur], ctx=cv.getContext("2d"), dpr=window.devicePixelRatio||1, w=cv.clientWidth, h=cv.clientHeight;
  if(cv.width!==Math.round(w*dpr)||cv.height!==Math.round(h*dpr)){ cv.width=Math.round(w*dpr); cv.height=Math.round(h*dpr); }
  ctx.setTransform(dpr,0,0,dpr,0,0); ctx.fillStyle="#ffffff"; ctx.fillRect(0,0,w,h);
  if(!v.topo){ ctx.fillStyle="#52514e"; ctx.font="13px system-ui,sans-serif"; ctx.fillText("No topology layout is stored for this version.",20,30); $("status").textContent=v.name; return; }
  const T=v.topo, G=topoGeometry(), fk=cur+"|"+$("tmode").value; if(TV.fitted!==fk){ topoLoad(); topoFit(); TV.fitted=fk; recompute=true; }
  if(recompute||!topoCache||topoCache.key!==fk){ const r=compute(); lastOut=r; const col=topoColours(r), key=$("tnode").value, meta=scoreMeta(), sc=v.scores[$("score").value]; let nodes=null, nbar=null;
    if(key!=="structure"){ const all=[]; nodes=T.isl.map((t,k)=>{ const np=nodeProperty(v,v.S[t.g],G.rows[k],r,sc,key,meta); np.val.forEach(x=>{ if(Number.isFinite(x)) all.push(x); }); nbar=np.title; return np.val; });
      const rg=ramp(all); nodes=nodes.map(a=>a.map(x=>rg.bin(x))); nbar={title:"nodes: "+nbar, lo:rg.lo, hi:rg.hi}; }
    const prow=G.rows.map(a=>Int32Array.from(a).sort((x,y)=>col.ci[x]-col.ci[y]));
    topoCache={key:fk,r,col,nodes,nbar,prow};
    $("outinfo").innerHTML = `<b>${r.nout.toLocaleString()}</b> outliers among ${r.nshown.toLocaleString()} shown points (${r.nshown? (100*r.nout/r.nshown).toFixed(2):"0"} %)<br><span class="kind">${meta.low?"low":"high"} values are outliers · ${r.info||""}</span>`; }
  const {r,col,nodes,nbar,prow}=topoCache, S_=TV.scale, K=T.isl.length;
  // affine map local -> screen for every island: screen = (A local + c) scaled; A = view rotation x island rotation
  const M=new Array(K), vis=new Array(K), ctr=new Array(K), PS=topoPersp(G), cs_=new Array(K);
  for(let k=0;k<K;k++){ vis[k]=topoVisible(k); const A=mm(TV.R,T.isl[k].R), c=mv(TV.R,topoCentre(k,G)); ctr[k]=c; cs_[k]=topoScreen(c,PS);
    M[k]=[S_*A[0],S_*A[1],S_*A[2],TV.tx+S_*c[0], -S_*A[3],-S_*A[4],-S_*A[5],TV.ty-S_*c[1], A[6],A[7],A[8],c[2]]; }
  // screen x of a local point; PY then holds its screen y and PF its perspective size factor
  let PF=1, PY=0;
  const sx=(m,L,i)=>{ const x=m[0]*L[3*i]+m[1]*L[3*i+1]+m[2]*L[3*i+2]+m[3], y=m[4]*L[3*i]+m[5]*L[3*i+1]+m[6]*L[3*i+2]+m[7]; if(!PS.on){ PF=1; PY=y; return x; }
    PF=PS.f(m[8]*L[3*i]+m[9]*L[3*i+1]+m[10]*L[3*i+2]+m[11]); PY=PS.cy+(y-PS.cy)*PF; return PS.cx+(x-PS.cx)*PF; };
  const sy=()=>PY;                                                   // always called right after sx for the same point
  const pt=(k,q)=>{ const x=sx(M[k],q,0); return [x,PY]; };
  const showL=$("tlabels").checked, hov=topoHover, showP=$("tpts").checked, only=$("only").checked, hl=$("hl").checked;
  const order=[...Array(K).keys()].filter(k=>vis[k]).sort((a,b)=>ctr[a][2]-ctr[b][2]);       // far structures first
  const psz=Math.max(1.2,Math.min(3,+$("psize").value*0.8)), pop=Math.min(1,+$("popac").value+0.1), pl=G.pl;
  order.forEach(k=>{ const g=v.S[T.isl[k].g], m=M[k], L=G.loc[k], c0=cs_[k], kf=c0[2];
    if(showL){ ctx.beginPath(); ctx.arc(c0[0],c0[1],G.R[k]*1.05*S_*kf,0,2*Math.PI); ctx.globalAlpha = k===hov ? 0.12 : 0.04; ctx.fillStyle=g.colour; ctx.fill();
      ctx.globalAlpha = k===hov ? 0.9 : 0.3; ctx.strokeStyle=g.colour; ctx.lineWidth = k===hov ? 1.6 : 0.8; ctx.stroke(); }
    if(showP && !only){ ctx.globalAlpha=pop; let last=-1; const rows=prow[k]; for(let q=0;q<rows.length;q++){ const i=rows[q]; if(!r.shown[i]||(hl&&r.out[i])) continue;
        if(col.ci[i]!==last){ last=col.ci[i]; ctx.fillStyle=col.pal[last]; } const x=sx(m,pl,i), z=psz*PF; ctx.fillRect(x-z/2, PY-z/2, z, z); } }
    const E=g.edges; ctx.strokeStyle=g.colour; ctx.globalAlpha=0.75; ctx.lineWidth=0.8; ctx.beginPath();
    for(let e=0;e<E.length/2;e++){ const a=E[2*e], b=E[2*e+1]; ctx.moveTo(sx(m,L,a),sy(m,L,a)); ctx.lineTo(sx(m,L,b),sy(m,L,b)); }
    ctx.stroke(); ctx.globalAlpha=1; const nr=Math.max(1.6,Math.min(4.5,0.022*G.R[k]*S_+1.4));
    for(let i=0;i<g.n_nodes;i++){ ctx.beginPath(); const x=sx(m,L,i); ctx.arc(x,PY,nr*PF,0,2*Math.PI);
      ctx.fillStyle = nodes ? (nodes[k][i]<0 ? "#e4e3dd" : RAMP[nodes[k][i]]) : g.colour; ctx.fill(); if(nodes){ ctx.strokeStyle="#ffffff"; ctx.lineWidth=0.5; ctx.stroke(); } }
    if(showP && (hl||only)){ ctx.fillStyle="#e34948"; ctx.strokeStyle="#ffffff"; ctx.lineWidth=0.5; const rows=G.rows[k];
      for(let q=0;q<rows.length;q++){ const i=rows[q]; if(!r.out[i]) continue; const x=sx(m,pl,i), y=PY, q4=4*PF; ctx.beginPath(); ctx.moveTo(x,y-q4); ctx.lineTo(x+q4,y); ctx.lineTo(x,y+q4); ctx.lineTo(x-q4,y); ctx.closePath(); ctx.fill(); ctx.stroke(); } }
  });
  if($("tlinks").checked){ const minn=+$("tlink").value; ctx.lineCap="round"; G.links.forEach(l=>{ if(l.n<minn||!vis[l.a]||!vis[l.b]) return; const A=pt(l.a,l.pa), B=pt(l.b,l.pb), mine = hov===l.a||hov===l.b;
      ctx.globalAlpha = hov<0 ? 0.4 : (mine ? 0.75 : 0.08); ctx.strokeStyle = mine ? "#0b0b0b" : "#52514e"; ctx.lineWidth=0.8+9*l.n/G.maxn; ctx.beginPath(); ctx.moveTo(A[0],A[1]); ctx.lineTo(B[0],B[1]); ctx.stroke(); }); ctx.globalAlpha=1; }
  if(showL){ ctx.font="12px system-ui,sans-serif"; ctx.textAlign="center"; ctx.textBaseline="bottom"; order.forEach(k=>{ const g=v.S[T.isl[k].g], t=g.label+(g.kind==="group"?" (rest)":"")+" · "+g.dim+"-D";
      const y=cs_[k][1]-G.R[k]*1.05*S_*cs_[k][2]-3, x=cs_[k][0], tw=ctx.measureText(t).width; ctx.fillStyle="rgba(255,255,255,0.8)"; ctx.fillRect(x-tw/2-3,y-14,tw+6,15); ctx.fillStyle="#0b0b0b"; ctx.fillText(t,x,y); }); }
  // orientation of the frame
  { const o=[44,h-44], len=26, cols=["#e34948","#1baf7a","#2a78d6"]; ctx.font="10px system-ui,sans-serif"; ctx.textAlign="center"; ctx.textBaseline="middle";
    [0,1,2].map(q=>({q,z:TV.R[6+q]})).sort((a,b)=>a.z-b.z).forEach(({q})=>{ const ex=TV.R[q], ey=TV.R[3+q]; ctx.strokeStyle=cols[q]; ctx.lineWidth=2; ctx.beginPath(); ctx.moveTo(o[0],o[1]); ctx.lineTo(o[0]+len*ex,o[1]-len*ey); ctx.stroke();
      ctx.fillStyle=cols[q]; ctx.fillText(String(q+1),o[0]+(len+8)*ex,o[1]-(len+8)*ey); }); }
  let by=14; [nbar,col.bar].forEach(b=>{ if(!b) return; const bw=16, x0=w-18-7*bw; ctx.textAlign="left"; ctx.textBaseline="top"; ctx.font="11px system-ui,sans-serif"; const tw=Math.max(ctx.measureText(b.title).width,7*bw);
    ctx.fillStyle="rgba(255,255,255,0.85)"; ctx.fillRect(w-26-tw,by-4,tw+16,46); ctx.fillStyle="#52514e"; ctx.fillText(b.title,w-18-tw,by); RAMP.forEach((c,i)=>{ ctx.fillStyle=c; ctx.fillRect(x0+i*bw,by+15,bw,9); });
    ctx.fillStyle="#52514e"; ctx.fillText(b.lo.toPrecision(3),x0,by+26); ctx.textAlign="right"; ctx.fillText(b.hi.toPrecision(3),x0+7*bw,by+26); by+=52; });
  const nv=order.length, nl=G.links.filter(l=>l.n>=+$("tlink").value&&vis[l.a]&&vis[l.b]).length;
  $("status").textContent = v.name+" · "+nv+" structures, "+nl+" links"+(hov>=0 ? " · structure "+v.S[T.isl[hov].g].label+": "+v.S[T.isl[hov].g].dim+"-D, "+v.S[T.isl[hov].g].n_nodes+" nodes, "+G.rows[hov].length.toLocaleString()+" points · drag to move, Shift-drag to turn" : (G.key.endsWith("flat") && T.gap ? " · mean gap between nodes in contact: "+T.gap.flat.toFixed(2)+" edge lengths (islands at joint positions "+T.gap.islands.toFixed(2)+", joint "+T.gap.joint.toFixed(2)+")" : " · drag the background to rotate the frame, drag a structure to move it"));
}
function topoRedraw(){ if(topoRAF) return; topoRAF=requestAnimationFrame(()=>{ topoRAF=0; if(view==="topo") drawTopo(false); }); }
(function(){ const cv=$("topo"), pos=e=>{ const b=cv.getBoundingClientRect(); return [e.clientX-b.left, e.clientY-b.top]; };
  cv.addEventListener("contextmenu",e=>e.preventDefault());
  cv.addEventListener("pointerdown",e=>{ if(!V[cur].topo) return; const [mx,my]=pos(e), pan=e.button===2||e.ctrlKey||e.metaKey, k=pan ? -1 : topoHit(mx,my); cv.setPointerCapture(e.pointerId);
    if(k>=0){ const t=V[cur].topo.isl[k]; const G=topoGeometry(); topoDrag={kind:e.shiftKey?"turn":"move", k, mx, my, d:t.d.slice(), R:t.R.slice(), f:topoPersp(G).f(mv(TV.R,topoCentre(k,G))[2])}; }
    else topoDrag={kind:pan?"pan":"orbit", mx, my, tx:TV.tx, ty:TV.ty, R:TV.R.slice()}; cv.style.cursor="grabbing"; e.preventDefault(); });
  cv.addEventListener("pointermove",e=>{ if(!V[cur].topo) return; const [mx,my]=pos(e);
    if(!topoDrag){ const k=topoHit(mx,my); if(k!==topoHover){ topoHover=k; cv.style.cursor = k>=0 ? "move" : "grab"; topoRedraw(); } return; }
    const dx=mx-topoDrag.mx, dy=my-topoDrag.my, D_=topoDrag;
    if(D_.kind==="pan"){ TV.tx=D_.tx+dx; TV.ty=D_.ty+dy; }
    else if(D_.kind==="orbit"){ // rotate the whole frame about the centre of the window, structures stay rigidly arranged
      const cx=cv.clientWidth/2, cy=cv.clientHeight/2, Rn=orthonormal(mm(rotXY(dy*0.008,dx*0.008),D_.R));
      const pivot=mv(mt(D_.R),[(cx-D_.tx)/TV.scale, -(cy-D_.ty)/TV.scale, 0]), q=mv(Rn,pivot); TV.R=Rn; TV.tx=cx-TV.scale*q[0]; TV.ty=cy+TV.scale*q[1]; }
    else { const t=V[cur].topo.isl[D_.k];
      if(D_.kind==="move"){ const wd=mv(mt(TV.R),[dx/(TV.scale*D_.f),-dy/(TV.scale*D_.f),0]); t.d=[D_.d[0]+wd[0],D_.d[1]+wd[1],D_.d[2]+wd[2]]; }
      else t.R=orthonormal(mm(mm(mt(TV.R),mm(rotXY(dy*0.01,dx*0.01),TV.R)),D_.R)); }
    topoRedraw(); });
  const end=()=>{ if(topoDrag && topoDrag.kind!=="pan") topoSave(); topoDrag=null; cv.style.cursor = topoHover>=0 ? "move" : "grab"; };
  cv.addEventListener("pointerup",end); cv.addEventListener("pointercancel",end);
  cv.addEventListener("pointerleave",()=>{ if(!topoDrag && topoHover>=0){ topoHover=-1; topoRedraw(); } });
  cv.addEventListener("wheel",e=>{ e.preventDefault(); const [mx,my]=pos(e), f=Math.exp(-e.deltaY*0.0015); TV.tx=mx-(mx-TV.tx)*f; TV.ty=my-(my-TV.ty)*f; TV.scale*=f; topoRedraw(); },{passive:false});
  cv.addEventListener("dblclick",e=>{ const [mx,my]=pos(e), k=topoHit(mx,my); if(k<0) return; const g=String(V[cur].topo.isl[k].g), sel=$("gstruct"); let o=[...sel.options].find(o=>o.value===g) || [...sel.options].find(o=>V[cur].S[+o.value] && String(V[cur].S[+o.value].graph_of)===g); if(o){ sel.value=o.value; setView("graph"); } });
  window.addEventListener("keydown",e=>{ if(view!=="topo"||topoHover<0||(e.key!=="f"&&e.key!=="F")||e.target.tagName==="INPUT"||e.target.tagName==="SELECT") return; const t=V[cur].topo.isl[topoHover]; t.R=mm(t.R,[-1,0,0,0,1,0,0,0,1]); topoSave(); topoRedraw(); });
  window.addEventListener("resize",()=>{ if(view==="topo") topoRedraw(); });
  $("tmode").addEventListener("input",()=>{ $("tsep").value = $("tmode").value==="isl" ? 1 : 0; TV.fitted=""; draw(); });
  $("tsep").addEventListener("input",topoRedraw); $("tsep").addEventListener("change",()=>{ topoFit(); topoRedraw(); });
  $("tlink").addEventListener("input",()=>{ $("tlinkv").textContent=$("tlink").value; topoRedraw(); });
  ["tpts","tlinks","tlabels","tpersp"].forEach(id=>$(id).addEventListener("input",topoRedraw)); $("tnode").addEventListener("input",()=>{ if(view==="topo") drawTopo(true); });
  $("treset").onclick=()=>{ V[cur].topo && V[cur].topo.isl.forEach(t=>{ t.d=[0,0,0]; t.R=I3(); }); TV.R=rotXY(-0.35,0.45); topoSave(); topoFit(); topoRedraw(); };
  $("tfit").onclick=()=>{ topoFit(); topoRedraw(); };
  $("tpng").onclick=()=>{ const a=document.createElement("a"); a.href=cv.toDataURL("image/png"); a.download="topology.png"; a.click(); };
})();
function setView(w){ view=w; document.body.className = w==="graph" ? "vgraph" : w==="topo" ? "vtopo" : "v3d"; document.querySelectorAll(".tab").forEach(t=>t.classList.toggle("active", t.dataset.view===w)); Plotly.purge("plot"); lastThree=null; draw(); }
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
  const only=$("only").checked, hl=$("hl").checked, sc=v.scores[$("score").value], xe=extraOf(cm), isIdx = cm==="idxh"||cm==="idxs";
  const P=[[],[],[]], O=[[],[],[]], pc=[], cd=[], ocd=[];
  for(let i=0;i<N;i++){
    if(!r.shown[i]) continue; const b=i*D, isOut=r.out[i]===1;
    if(isOut && (hl||only)){ for(let a=0;a<nd;a++) O[a].push(X[b+ax[a]]); ocd.push([i,v.S[v.sid[i]].label,sc[i]]); if(hl) continue; }
    if(only) continue;
    for(let a=0;a<nd;a++) P[a].push(X[b+ax[a]]); cd.push([i,v.S[v.sid[i]].label,sc[i]]);
    if(cm==="structure") pc.push(v.S[v.sid[i]].colour); else if(cm==="dim") pc.push(DIMCOL[Math.min(v.S[v.sid[i]].dim,8)]); else if(isIdx) pc.push(v[cm] ? DIMCOL[Math.min(v[cm][i],8)] : DIMCOL[0]);
    else if(xe) pc.push(xe.kind==="num" ? xe.values[i] : (xe.values[i]<0 ? "#b9b8b2" : colourOf(xe.values[i]))); else pc.push(v.scores[cm][i]);
  }
  const traces=[], put=(t,A)=>{ t.x=A[0]; t.y=A[1]; if(three) t.z=A[2]; return t; };
  const ht="row %{customdata[0]} · %{customdata[1]}<br>score %{customdata[2]:.4g}<extra></extra>";
  const dens=densityTraces(ax);
  if(!three) dens.traces.forEach(t=>traces.push(t));
  if($("showpts").checked && P[0].length){
    const marker={size:+$("psize").value+(three?0:1), opacity:+$("popac").value, color:pc, line:{width:0}};
    if(cm!=="structure" && cm!=="dim" && !isIdx && !(xe && xe.kind==="cat")){ const rg=seqColour(pc); Object.assign(marker,{colorscale:"Blues",reversescale:true,cmin:rg.cmin,cmax:rg.cmax,
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
function downloadTable(){
  const v=V[cur], r=lastOut||compute(), keys=[...document.querySelectorAll("#tabcrit input:checked")].map(e=>e.value);
  const nan=new Set([...document.querySelectorAll("#tabnan input:checked")].map(e=>+e.dataset.code));
  const ws=$("tabstruct").checked, wsel=$("tabsel").checked, wx=$("tabx").checked, we=$("tabextra").checked, q=t=>'"'+String(t).replace(/"/g,'""')+'"';
  const head=["row","status"]; if(ws) head.push("structure","structure_dim","node","dim_index_raw","dim_index_smoothed"); if(wsel) head.push("shown","selected_outlier");
  keys.forEach(k=>head.push(k)); if(wx) DATA.names.forEach(t=>head.push(q(t))); if(we) EXTRA.forEach(e=>head.push(q(e.name)));
  const rows=[head.join(",")];
  for(let i=0;i<N;i++){ const st=v.status? v.status[i] : -1, blank=nan.has(st), s=v.S[v.sid[i]], c=[i, q(st<0? "unknown" : v.status_names[st])];
    if(ws) c.push(s? q(s.label):"", s? s.dim:"", v.node[i], v.idxh? v.idxh[i]:"", v.idxs? v.idxs[i]:"");
    if(wsel) c.push(r.shown[i]?1:0, r.out[i]?1:0);
    keys.forEach(k=>{ const t=v.scores[k][i]; c.push(blank||!Number.isFinite(t) ? "NaN" : t); });
    if(wx) for(let k=0;k<D;k++) c.push(X[i*D+k].toPrecision(7));
    if(we) EXTRA.forEach(e=>c.push(e.kind==="num" ? (Number.isFinite(e.values[i])? e.values[i] : "NaN") : (e.values[i]<0? "" : q(e.cats[e.values[i]]))));
    rows.push(c.join(",")); }
  const a=document.createElement("a"); a.href=URL.createObjectURL(new Blob([rows.join("\n")],{type:"text/csv"})); a.download="criteria_all_points.csv"; a.click(); URL.revokeObjectURL(a.href);
}
function setAll(fn){ const v=V[cur]; v.S.forEach(s=>{ s.on=fn(s,v); }); syncParents(v); buildTree(); draw(); }
$("all").onclick=()=>setAll(()=>true); $("none").onclick=()=>setAll(()=>false);
$("leavesonly").onclick=()=>setAll(s=>s.kind==="leaf");
$("top").onclick=()=>setAll(s=>s.kind==="residual");
$("version").onchange=()=>{ cur=+$("version").value; buildTree(); fillGraphSelect(); setThrRange(); draw(); };
document.querySelectorAll(".tab").forEach(t=>t.addEventListener("click",()=>setView(t.dataset.view)));
refreshExtraOptions();
["gstruct","gcolour","glayout","gpts","gfill","gloc","gy","gedge","gproj","gsize","isomode","isores"].forEach(id=>$(id).addEventListener("input",draw));
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
DATA.scores.forEach(s=>{ const b=document.createElement("b"); b.textContent=s.label+" ("+(s.low?"low":"high")+" values are outliers)"; $("scorehelp").append(b, document.createTextNode(s.help||"")); });
DATA.scores.forEach(s=>{ const l=document.createElement("label"); l.className="check"; l.title=s.help||""; const c=document.createElement("input"); c.type="checkbox"; c.value=s.key; c.checked=true; l.append(c," "+s.label+" ["+s.key+"]"); $("tabcrit").appendChild(l); });
$("tabdl").onclick=downloadTable;
$("csv").onclick=()=>download(true); $("csvall").onclick=()=>download(false);
$("lochead").addEventListener("click",()=>{ $("loc").classList.toggle("big"); Plotly.Plots.resize("locplot"); });
document.body.className="v3d"; buildTree(); fillGraphSelect(); setThrRange(); draw();
</script></body></html>
"""
