"""Static matplotlib figures for every stage of the pipeline.

All functions return a ``matplotlib.figure.Figure``. Data with more than three features are
shown through ``axes=(i, j, k)`` (feature indices, two or three of them) or, with
``projection="pca"``, in the leading principal components of the plotted points.
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "View", "plot_data", "plot_preprocess", "plot_density_counts", "plot_dim_index", "plot_simplex",
    "plot_confusion", "plot_graphs", "plot_crawl_progress", "plot_unfolded", "plot_agtm", "plot_loglik",
    "plot_curvature", "plot_node_values", "plot_filaments", "plot_profile", "filament_table",
    "plot_cross_sections", "plot_summary", "CATEGORICAL", "DIM_COLORS",
]

# validated categorical palette (fixed order) and neutral inks
CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
DIM_COLORS = {0: "#b9b8b2", 1: "#e34948", 2: "#2a78d6", 3: "#1baf7a", 4: "#eda100", 5: "#6250d6"}
INK, MUTED, GRID, NOISE = "#0b0b0b", "#52514e", "#dddcd6", "#b9b8b2"
SEQ = "Blues"


def _plt():
    import matplotlib.pyplot as plt

    return plt


def color_of(k: int) -> str:
    return CATEGORICAL[k % len(CATEGORICAL)]


class View:
    """Maps d-dimensional points to the 2 or 3 plotted coordinates."""

    def __init__(self, X, axes=None, projection: str | None = None, names=None):
        X = np.asarray(X, float)
        d = X.shape[1]
        self.names = list(names) if names is not None else [f"x{k + 1}" for k in range(d)]
        self.mean, self.basis = None, None
        if projection == "pca" and d > 2:
            self.mean = X.mean(0)
            _, _, vt = np.linalg.svd(X[np.random.default_rng(0).choice(len(X), min(len(X), 20000), replace=False)]
                                     - self.mean, full_matrices=False)
            k = 3 if axes is None else len(axes)
            self.basis = vt[:min(k, d)]
            self.labels = [f"PC{i + 1}" for i in range(len(self.basis))]
        else:
            self.axes = tuple(axes) if axes is not None else tuple(range(min(d, 3)))
            if not 2 <= len(self.axes) <= 3 or max(self.axes) >= d:
                raise ValueError("axes must list two or three valid feature indices")
            self.labels = [self.names[a] for a in self.axes]
        self.ndim = len(self.labels)

    def __call__(self, X):
        X = np.asarray(X, float)
        if self.basis is not None:
            return (X - self.mean) @ self.basis.T
        return X[:, self.axes]

    def new_axes(self, fig, *pos):
        ax = fig.add_subplot(*pos, projection="3d" if self.ndim == 3 else None)
        self.style(ax)
        return ax

    def style(self, ax):
        ax.set_xlabel(self.labels[0], color=MUTED)
        ax.set_ylabel(self.labels[1], color=MUTED)
        if self.ndim == 3:
            ax.set_zlabel(self.labels[2], color=MUTED)
            for a in (ax.xaxis, ax.yaxis, ax.zaxis):
                a.pane.set_facecolor("white")
                a.pane.set_edgecolor(GRID)
                a._axinfo["grid"]["color"] = GRID
                a._axinfo["grid"]["linewidth"] = 0.5
            try:
                ax.set_box_aspect(None)
            except Exception:
                pass
        else:
            ax.grid(color=GRID, linewidth=0.5)
            for s in ax.spines.values():
                s.set_color(GRID)
        ax.tick_params(colors=MUTED, labelsize=7)

    def equal(self, ax, P, pad=1.02, robust=None):
        """Equal aspect and limits from the plotted coordinates ``P`` (``robust``: quantile
        range used instead of the full extent, e.g. 0.005)."""
        if not len(P):
            return
        if robust:
            lo, hi = np.quantile(P, robust, axis=0), np.quantile(P, 1 - robust, axis=0)
        else:
            lo, hi = P.min(0), P.max(0)
        c, h = (lo + hi) / 2, (hi - lo).max() / 2 * pad + 1e-12
        if self.ndim == 2:
            ax.set_xlim(c[0] - h, c[0] + h)
            ax.set_ylim(c[1] - h, c[1] + h)
            ax.set_aspect("equal", adjustable="box")
        if self.ndim == 3:
            ax.set_xlim(c[0] - h, c[0] + h)
            ax.set_ylim(c[1] - h, c[1] + h)
            ax.set_zlim(c[2] - h, c[2] + h)
            ax.set_box_aspect((1, 1, 1))


def _view(X, axes, projection, names, view=None) -> View:
    return view if view is not None else View(X, axes, projection, names)


def _sub(n, max_points, seed=0):
    if max_points is None or n <= max_points:
        return slice(None)
    return np.sort(np.random.default_rng(seed).choice(n, max_points, replace=False))


def _scatter(ax, P, s=1.0, **kw):
    kw.setdefault("linewidths", 0)
    kw.setdefault("rasterized", True)
    return ax.scatter(*P.T, s=s, **kw)


def _edges(ax, view, graph, color=INK, lw=0.8, values=None, cmap=SEQ, vmin=None, vmax=None, alpha=1.0):
    P = view(graph.nodes)
    seg = np.stack([P[graph.edges[:, 0]], P[graph.edges[:, 1]]], 1)
    if view.ndim == 3:
        from mpl_toolkits.mplot3d.art3d import Line3DCollection as LC
    else:
        from matplotlib.collections import LineCollection as LC
    if values is None:
        lc = LC(seg, colors=color, linewidths=lw, alpha=alpha)
    else:
        lc = LC(seg, cmap=cmap, linewidths=lw, alpha=alpha)
        lc.set_array(np.asarray(values, float))
        lc.set_clim(vmin, vmax)
    ax.add_collection(lc) if view.ndim == 2 else ax.add_collection3d(lc)
    if view.ndim == 2:
        ax.update_datalim(P)
        ax.autoscale_view()
    return lc


def _title(ax, text):
    ax.set_title(text, fontsize=9, color=INK, loc="left")


def _legend(ax, entries, **kw):
    from matplotlib.lines import Line2D

    h = [Line2D([], [], marker="o", ls="", color=c, markersize=5, label=l) for l, c in entries]
    kw.setdefault("fontsize", 7)
    kw.setdefault("frameon", False)
    return ax.legend(handles=h, **kw)


# ------------------------------------------------------------------ data & preprocessing
def plot_data(X, color=None, *, labels=None, axes=None, projection=None, names=None, max_points=60000,
              title="Data", color_label=None, figsize=(6, 5)):
    """Point cloud, optionally coloured by a continuous value (``color``) or by integer ``labels``
    (0 is drawn as background)."""
    plt = _plt()
    X = np.asarray(X, float)
    v = View(X, axes, projection, names)
    fig = plt.figure(figsize=figsize)
    ax = v.new_axes(fig, 1, 1, 1)
    sel = _sub(len(X), max_points)
    P = v(X[sel])
    if labels is not None:
        lab = np.asarray(labels)[sel]
        _scatter(ax, P[lab == 0], 0.4, color=NOISE, alpha=0.5)
        ents = []
        for i, k in enumerate(np.unique(lab[lab != 0])):
            _scatter(ax, P[lab == k], 1.2, color=color_of(i))
            ents.append((str(k), color_of(i)))
        if 1 < len(ents) <= 12:
            _legend(ax, ents, loc="upper right")
    elif color is not None:
        sc = _scatter(ax, P, 1.2, c=np.asarray(color, float)[sel], cmap=SEQ)
        cb = fig.colorbar(sc, ax=ax, shrink=0.7, pad=0.08)
        cb.set_label(color_label or "", fontsize=8, color=MUTED)
    else:
        _scatter(ax, P, 0.8, color=CATEGORICAL[0], alpha=0.6)
    v.equal(ax, P)
    _title(ax, f"{title} ({len(X):,} points)")
    fig.tight_layout()
    return fig


def plot_preprocess(X, result, *, axes=None, projection=None, names=None, max_points=60000, figsize=(13, 4.2)):
    """Raw data with removed points greyed out, the diffused points, and how far points moved."""
    plt = _plt()
    X = np.asarray(X, float)
    v = View(X, axes, projection, names)
    fig = plt.figure(figsize=figsize)
    kept = np.zeros(len(X), bool)
    kept[result.indices] = True
    ax = v.new_axes(fig, 1, 3, 1)
    s = _sub(len(X), max_points)
    P, k = v(X[s]), kept[s]
    _scatter(ax, P[~k], 0.4, color=NOISE, alpha=0.5)
    _scatter(ax, P[k], 0.8, color=CATEGORICAL[0])
    v.equal(ax, P)
    _title(ax, f"Noisy data: kept {kept.sum():,} of {len(X):,}")
    _legend(ax, [("kept", CATEGORICAL[0]), ("removed", NOISE)], loc="upper right")
    ax2 = v.new_axes(fig, 1, 3, 2)
    s2 = _sub(len(result), max_points)
    P2 = v(result.diffused[s2])
    _scatter(ax2, P2, 0.8, color=CATEGORICAL[1])
    v.equal(ax2, P)
    _title(ax2, "Diffused points")
    ax3 = fig.add_subplot(1, 3, 3)
    disp = np.linalg.norm(result.diffused - result.noisy, axis=1)
    ax3.hist(disp, bins=60, color=CATEGORICAL[0], edgecolor="white", linewidth=0.3)
    ax3.set_xlabel("displacement of a point", color=MUTED)
    ax3.set_ylabel("points", color=MUTED)
    ax3.grid(axis="y", color=GRID, linewidth=0.5)
    ax3.set_axisbelow(True)
    for sp in ("top", "right"):
        ax3.spines[sp].set_visible(False)
    ax3.tick_params(colors=MUTED, labelsize=7)
    _title(ax3, f"Displacement (median {np.median(disp):.3g})")
    fig.tight_layout()
    return fig


def plot_density_counts(X, radius, min_points=None, figsize=(5.5, 3.6)):
    """Histogram of neighbour counts within ``radius``; helps choosing the density threshold."""
    from ._util import as_points, build_tree

    plt = _plt()
    X = as_points(X)
    c = build_tree(X).radius_count(X, float(radius))
    fig, ax = plt.subplots(figsize=figsize)
    ax.hist(c, bins=np.arange(0.5, min(c.max(), np.quantile(c, 0.995) + 1) + 1.5), color=CATEGORICAL[0],
            edgecolor="white", linewidth=0.2)
    if min_points is not None:
        ax.axvline(min_points - 0.5, color=INK, linewidth=1, linestyle="--")
        ax.annotate(f"threshold = {min_points}", (min_points - 0.5, ax.get_ylim()[1] * 0.92), fontsize=7,
                    color=INK, xytext=(4, 0), textcoords="offset points")
    ax.set_xlabel(f"points within radius {radius:g} (including the point itself)", color=MUTED)
    ax.set_ylabel("points", color=MUTED)
    ax.grid(axis="y", color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(colors=MUTED, labelsize=7)
    _title(ax, "Neighbourhood sizes")
    fig.tight_layout()
    return fig


# ------------------------------------------------------------------ dimensionality index
def plot_dim_index(points, result, *, axes=None, projection=None, names=None, max_points=60000,
                   per_dimension=True, figsize=None):
    """Points coloured by dimensionality index: hard and smoothed side by side, and (optionally)
    one panel per dimension with its number of points."""
    plt = _plt()
    P_all = np.asarray(points, float)
    v = View(P_all, axes, projection, names)
    dims = [j for j in np.unique(result.index) if j > 0]
    rows = 2 if per_dimension and dims else 1
    ncol = max(2, len(dims))
    fig = plt.figure(figsize=figsize or (4.3 * ncol, 4.0 * rows))
    s = _sub(len(P_all), max_points)
    P = v(P_all[s])
    for c, (name, idx) in enumerate((("Hard index", result.index_hard), ("Smoothed index", result.index))):
        ax = v.new_axes(fig, rows, ncol, c + 1)
        idx = np.asarray(idx)[s]
        for j in np.unique(idx):
            _scatter(ax, P[idx == j], 0.8, color=DIM_COLORS.get(int(j), INK))
        v.equal(ax, P)
        _title(ax, name)
        _legend(ax, [("undetermined" if j == 0 else f"{j}-D", DIM_COLORS.get(int(j), INK))
                     for j in np.unique(result.index)], loc="upper right")
    if rows == 2:
        idx = np.asarray(result.index)[s]
        for c, j in enumerate(dims):
            ax = v.new_axes(fig, rows, ncol, ncol + c + 1)
            _scatter(ax, P[idx == j], 0.8, color=DIM_COLORS.get(int(j), INK))
            v.equal(ax, P)
            _title(ax, f"{j}-D: {(result.index == j).sum():,} points")
    fig.tight_layout()
    return fig


def plot_simplex(result, max_points=30000, figsize=(9, 4)):
    """Normalised local spectra and their barycentric coordinates (first three components),
    coloured by the smoothed index (paper Fig. 11)."""
    plt = _plt()
    fig = plt.figure(figsize=figsize)
    s = _sub(len(result.index), max_points)
    idx = result.index[s]
    for c, (A, lab, name) in enumerate(((result.eigenvalues, "λ", "Normalised spectra"),
                                         (result.barycentric, "α", "Barycentric coordinates"))):
        A = A[s]
        three = A.shape[1] >= 3
        ax = fig.add_subplot(1, 2, c + 1, projection="3d" if three else None)
        for j in np.unique(idx):
            _scatter(ax, A[idx == j][:, :3 if three else 2], 1.0, color=DIM_COLORS.get(int(j), INK))
        ax.set_xlabel(f"{lab}1", color=MUTED)
        ax.set_ylabel(f"{lab}2", color=MUTED)
        if three:
            ax.set_zlabel(f"{lab}3", color=MUTED)
            ax.view_init(25, 35)
        ax.tick_params(colors=MUTED, labelsize=7)
        _title(ax, name)
    fig.tight_layout()
    return fig


def plot_confusion(index, true_dim, figsize=(4.2, 3.8)):
    """Column-normalised confusion matrix of output index against true dimension (paper Fig. 12)."""
    plt = _plt()
    index, true_dim = np.asarray(index), np.asarray(true_dim)
    classes = np.union1d(np.unique(index), np.unique(true_dim))
    M = np.array([[np.sum((index == a) & (true_dim == b)) for b in classes] for a in classes], float)
    M /= np.maximum(M.sum(0, keepdims=True), 1)
    fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(M, cmap="Greys", vmin=0, vmax=1)
    name = lambda k: "noise / none" if k == 0 else f"D = {k}"  # noqa: E731
    ax.set_xticks(range(len(classes)), [name(k) for k in classes], fontsize=7)
    ax.set_yticks(range(len(classes)), [name(k) for k in classes], fontsize=7)
    for i in range(len(classes)):
        for j in range(len(classes)):
            ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=8,
                    color="white" if M[i, j] > 0.5 else INK)
    ax.set_xlabel("true", color=MUTED)
    ax.set_ylabel("estimated", color=MUTED)
    _title(ax, "Dimensionality index")
    fig.tight_layout()
    return fig


# ------------------------------------------------------------------ graphs
def _graph_list(obj):
    """Accept a learner, a dict dim -> graphs, a list of graphs / Manifold objects."""
    if hasattr(obj, "manifolds_"):
        return [m.graph for m in obj.manifolds_]
    if isinstance(obj, dict):
        return [g for gs in obj.values() for g in gs]
    return [getattr(g, "graph", g) for g in obj]


def plot_graphs(graphs, data=None, *, axes=None, projection=None, names=None, max_points=40000,
                panels=False, number=True, figsize=None):
    """Embedded graphs, one colour per manifold, over the (greyed) data. ``panels=True`` adds one
    panel per manifold."""
    plt = _plt()
    graphs = _graph_list(graphs)
    ref = np.asarray(data, float) if data is not None else np.vstack([g.nodes for g in graphs])
    v = View(ref, axes, projection, names)
    n = len(graphs)
    ncol = 1 if not panels else min(4, n + 1)
    nrow = 1 if not panels else int(np.ceil((n + 1) / ncol))
    fig = plt.figure(figsize=figsize or ((6.5, 5.5) if not panels else (4 * ncol, 3.8 * nrow)))
    ax = v.new_axes(fig, nrow, ncol, 1)
    P = v(ref[_sub(len(ref), max_points)])
    if data is not None:
        _scatter(ax, P, 0.3, color=NOISE, alpha=0.4)
    for k, g in enumerate(graphs):
        _edges(ax, v, g, color_of(k), 0.9)
        if number and n > 1:
            c = v(g.nodes).mean(0)
            ax.text(*c, str(k), fontsize=8, color=INK, ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec=GRID, lw=0.5))
    v.equal(ax, v(np.vstack([g.nodes for g in graphs])) if graphs else P, pad=1.12)
    _title(ax, f"{n} manifold graph{'s' if n != 1 else ''}")
    if panels:
        for k, g in enumerate(graphs):
            a = v.new_axes(fig, nrow, ncol, k + 2)
            if data is not None and g.points is not None and len(g.points) <= len(ref):
                pass
            _edges(a, v, g, color_of(k), 0.9)
            Q = v(g.nodes)
            _scatter(a, Q, 4, color=color_of(k))
            v.equal(a, Q)
            _title(a, f"{k}: {g.dim}-D, {g.n_nodes} nodes")
    fig.tight_layout()
    return fig


def plot_crawl_progress(X, r, dim, *, seed=0, steps=(5, 20, 80, 320, None), eta=0.75, beta=0.4,
                        axes=None, projection=None, names=None, figsize=None):
    """Snapshots of the crawl after a growing number of nodes (paper Fig. 8)."""
    from .crawling import crawl

    plt = _plt()
    X = np.asarray(X, float)
    v = View(X, axes, projection, names)
    fig = plt.figure(figsize=figsize or (3.4 * len(steps), 3.4))
    P = v(X[_sub(len(X), 20000)])
    for c, mx in enumerate(steps):
        g = crawl(X, r, dim, seed=seed, eta=eta, beta=beta, max_nodes=mx)
        ax = v.new_axes(fig, 1, len(steps), c + 1)
        _scatter(ax, P, 0.3, color=NOISE, alpha=0.4)
        _edges(ax, v, g, CATEGORICAL[0], 0.8)
        v.equal(ax, P)
        _title(ax, f"{g.n_nodes} nodes")
    fig.tight_layout()
    return fig


def plot_unfolded(graph, values=None, *, node_values=None, label="", figsize=(5.5, 5), seed=0):
    """Planar layout of the abstract graph (paper Figs. 10 and 23), edges optionally coloured by
    ``values`` (one per edge, e.g. curvature) or nodes by ``node_values``."""
    import networkx as nx
    from matplotlib.collections import LineCollection

    plt = _plt()
    graph = getattr(graph, "graph", graph)
    G = nx.Graph()
    G.add_nodes_from(range(graph.n_nodes))
    G.add_edges_from(graph.edges.tolist())
    pos = nx.kamada_kawai_layout(G) if graph.n_nodes <= 1500 else nx.spring_layout(G, seed=seed, iterations=60)
    P = np.array([pos[k] for k in range(graph.n_nodes)])
    fig, ax = plt.subplots(figsize=figsize)
    seg = np.stack([P[graph.edges[:, 0]], P[graph.edges[:, 1]]], 1)
    if values is not None:
        lc = LineCollection(seg, cmap=SEQ, linewidths=1.6)
        lc.set_array(np.asarray(values, float))
        fig.colorbar(lc, ax=ax, shrink=0.7).set_label(label, fontsize=8, color=MUTED)
    else:
        lc = LineCollection(seg, colors=MUTED, linewidths=0.8)
    ax.add_collection(lc)
    if node_values is not None:
        sc = ax.scatter(*P.T, s=14, c=np.asarray(node_values, float), cmap=SEQ, zorder=3, edgecolors="white",
                        linewidths=0.4)
        fig.colorbar(sc, ax=ax, shrink=0.7).set_label(label, fontsize=8, color=MUTED)
    ax.update_datalim(P)
    ax.autoscale_view()
    ax.set_aspect("equal")
    ax.axis("off")
    _title(ax, f"Abstract graph ({graph.n_nodes} nodes, {graph.n_edges} edges)")
    fig.tight_layout()
    return fig


# ------------------------------------------------------------------ AGTM
def _models(obj):
    if hasattr(obj, "manifolds_"):
        return [m.agtm for m in obj.manifolds_ if m.agtm is not None]
    if hasattr(obj, "models"):
        return list(obj.models)
    if hasattr(obj, "means_"):
        return [obj]
    return [getattr(m, "agtm", m) for m in obj]


def _ellipse_outline(mean, cov, n_std, v, n=40):
    w, U = np.linalg.eigh(cov)
    out = []
    t = np.linspace(0, 2 * np.pi, n)
    circ = np.stack([np.cos(t), np.sin(t)], 1)
    d = len(mean)
    pairs = [(d - 1, d - 2)] if v.ndim == 2 else [(d - 1, d - 2), (d - 1, max(d - 3, 0)), (d - 2, max(d - 3, 0))]
    for a, b in pairs:
        if a == b:
            continue
        pts = mean + n_std * (circ[:, :1] * np.sqrt(max(w[a], 0)) * U[:, a] + circ[:, 1:] * np.sqrt(max(w[b], 0)) * U[:, b])
        out.append(v(pts))
    return out


def plot_agtm(models, data=None, *, axes=None, projection=None, names=None, max_points=40000,
              ellipsoids=120, n_std=1.5, figsize=(7, 6)):
    """Fitted AGTMs: embedded graphs (fitted node means) with a sample of covariance ellipsoids
    showing the manifold-aligned noise model, over the data."""
    plt = _plt()
    models = _models(models)
    ref = np.asarray(data, float) if data is not None else np.vstack([m.means_ for m in models])
    v = View(ref, axes, projection, names)
    fig = plt.figure(figsize=figsize)
    ax = v.new_axes(fig, 1, 1, 1)
    P = v(ref[_sub(len(ref), max_points)])
    if data is not None:
        _scatter(ax, P, 0.3, color=NOISE, alpha=0.35)
    total = sum(len(m.means_) for m in models)
    rng = np.random.default_rng(0)
    from .crawling import ManifoldGraph

    for k, m in enumerate(models):
        g = ManifoldGraph(m.means_, m.graph.edges, m.graph.dim)
        _edges(ax, v, g, color_of(k), 0.9)
        n_e = max(1, int(round(ellipsoids * len(m.means_) / total))) if ellipsoids else 0
        for i in rng.choice(len(m.means_), min(n_e, len(m.means_)), replace=False):
            for outline in _ellipse_outline(m.means_[i], m.covariances_[i], n_std, v):
                ax.plot(*outline.T, color=color_of(k), linewidth=0.5, alpha=0.7)
        bad = m.diagnostics_["degenerate"]
        if len(bad):
            _scatter(ax, v(m.means_[bad]), 18, color=INK, marker="x", linewidths=0.8, rasterized=False)
    v.equal(ax, v(np.vstack([m.means_ for m in models])), pad=1.12)
    _title(ax, f"AGTM: {len(models)} model{'s' if len(models) != 1 else ''}, {total} components "
               f"(ellipses: {n_std:g} std; x: degenerate)")
    fig.tight_layout()
    return fig


def plot_loglik(models, figsize=(5.5, 3.6)):
    """Mean log-likelihood per point along the EM iterations of each model."""
    plt = _plt()
    models = _models(models)
    fig, ax = plt.subplots(figsize=figsize)
    for k, m in enumerate(models):
        ax.plot(m.loglik_, color=color_of(k), linewidth=1.6, label=f"{k}: {m.graph.dim}-D, {len(m.means_)} nodes")
        for ev in m.diagnostics_["events"]:
            it = min(ev["iteration"], len(m.loglik_) - 1)
            ax.plot(it, m.loglik_[it], marker="x", color=INK, markersize=6)
    ax.set_xlabel("EM iteration", color=MUTED)
    ax.set_ylabel("mean log-likelihood per point", color=MUTED)
    ax.grid(color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(colors=MUTED, labelsize=7)
    if 1 < len(models) <= 10:
        ax.legend(fontsize=7, frameon=False)
    _title(ax, "EM convergence (x: degeneracy event)")
    fig.tight_layout()
    return fig


def plot_curvature(model, data=None, *, axes=None, projection=None, names=None, max_points=40000,
                   robust=0.98, figsize=(7, 6)):
    """Embedded graph with edges coloured by the directional curvature of the AGTM mapping."""
    plt = _plt()
    model = getattr(model, "agtm", model)
    from .crawling import ManifoldGraph

    g = ManifoldGraph(model.means_, model.graph.edges, model.graph.dim)
    c = model.edge_curvature()
    ref = np.asarray(data, float) if data is not None else g.nodes
    v = View(ref, axes, projection, names)
    fig = plt.figure(figsize=figsize)
    ax = v.new_axes(fig, 1, 1, 1)
    P = v(ref[_sub(len(ref), max_points)])
    if data is not None:
        _scatter(ax, P, 0.3, color=NOISE, alpha=0.35)
    lc = _edges(ax, v, g, values=c, lw=1.6, vmin=0, vmax=float(np.quantile(c, robust)))
    fig.colorbar(lc, ax=ax, shrink=0.7, pad=0.08).set_label("curvature", fontsize=8, color=MUTED)
    v.equal(ax, v(g.nodes), pad=1.12)
    _title(ax, f"Curvature along edges (median {np.median(c):.3g})")
    fig.tight_layout()
    return fig


def plot_node_values(model, X, values, *, label="", axes=None, projection=None, names=None, figsize=(7, 6)):
    """Graph nodes coloured by the responsibility-weighted mean of a per-point quantity (eq. 32)."""
    plt = _plt()
    model = getattr(model, "agtm", model)
    from .crawling import ManifoldGraph

    nv = model.node_means(X, values)
    g = ManifoldGraph(model.means_, model.graph.edges, model.graph.dim)
    v = View(np.asarray(X, float), axes, projection, names)
    fig = plt.figure(figsize=figsize)
    ax = v.new_axes(fig, 1, 1, 1)
    _edges(ax, v, g, MUTED, 0.6)
    Q = v(g.nodes)
    sc = _scatter(ax, Q, 22, c=nv, cmap=SEQ, edgecolors="white", linewidths=0.4, rasterized=False)
    fig.colorbar(sc, ax=ax, shrink=0.7, pad=0.08).set_label(label, fontsize=8, color=MUTED)
    v.equal(ax, Q)
    _title(ax, f"Node means of {label or 'values'}")
    fig.tight_layout()
    return fig


# ------------------------------------------------------------------ 1-D filaments
def chain_order(graph):
    """Order the nodes of a 1-D graph along the filament. Returns ``(order, arc_length)``; for
    branched graphs the longest path between leaves is used."""
    from scipy.sparse import csgraph

    graph = getattr(graph, "graph", graph)
    A = graph.adjacency(weighted=True)
    deg = graph.degree
    start = int(np.flatnonzero(deg <= 1)[0]) if (deg <= 1).any() else 0
    d0 = csgraph.dijkstra(A, directed=False, indices=start)
    a = int(np.nanargmax(np.where(np.isfinite(d0), d0, -1)))
    da, pred = csgraph.dijkstra(A, directed=False, indices=a, return_predecessors=True)
    b = int(np.nanargmax(np.where(np.isfinite(da), da, -1)))
    path = [b]
    while path[-1] != a and pred[path[-1]] >= 0:
        path.append(int(pred[path[-1]]))
    order = np.array(path[::-1])
    return order, da[order]


def filament_table(graphs) -> list[dict]:
    """Length, node count and number of attached points of each 1-D graph."""
    rows = []
    for k, g in enumerate(_graph_list(graphs)):
        if g.dim != 1:
            continue
        rows.append({"filament": k, "nodes": g.n_nodes, "length": g.total_length(),
                     "points": None if g.points is None else len(g.points),
                     "closed": bool((g.degree == 2).all())})
    return rows


def plot_filaments(graphs, data=None, *, axes=None, projection=None, names=None, max_points=40000,
                   max_panels=11, figsize=None):
    """All 1-D manifolds found, with distinct colours, plus one panel per filament showing its
    ordered node chain."""
    plt = _plt()
    fil = [g for g in _graph_list(graphs) if g.dim == 1]
    if not fil:
        raise ValueError("no 1-D manifolds to plot")
    ref = np.asarray(data, float) if data is not None else np.vstack([g.nodes for g in fil])
    v = View(ref, axes, projection, names)
    shown = fil[:max_panels]
    ncol = min(4, len(shown) + 1)
    nrow = int(np.ceil((len(shown) + 1) / ncol))
    fig = plt.figure(figsize=figsize or (4 * ncol, 3.8 * nrow))
    ax = v.new_axes(fig, nrow, ncol, 1)
    P = v(ref[_sub(len(ref), max_points)])
    if data is not None:
        _scatter(ax, P, 0.3, color=NOISE, alpha=0.4)
    for k, g in enumerate(fil):
        _edges(ax, v, g, color_of(k), 1.4)
    v.equal(ax, v(np.vstack([g.nodes for g in fil])), pad=1.25)
    _title(ax, f"{len(fil)} filament{'s' if len(fil) != 1 else ''} found")
    for k, g in enumerate(shown):
        a = v.new_axes(fig, nrow, ncol, k + 2)
        order, _ = chain_order(g)
        Q = v(g.nodes)
        _edges(a, v, g, color_of(k), 1.0)
        a.plot(*Q[order].T, color=INK, linewidth=0.8, marker="D", markersize=2.5)
        v.equal(a, Q)
        _title(a, f"Filament {k}: length {g.total_length():.3g}, {g.n_nodes} nodes")
    fig.tight_layout()
    return fig


def plot_profile(model, X=None, values=None, *, label="value", figsize=(6.5, 3.4)):
    """Profile along a filament: the responsibility-weighted mean of ``values`` at each node
    (or, without ``values``, the curvature) against arc length."""
    plt = _plt()
    model = getattr(model, "agtm", model)
    from .crawling import ManifoldGraph

    g = ManifoldGraph(model.means_, model.graph.edges, model.graph.dim)
    order, s = chain_order(g)
    if values is None:
        c = model.edge_curvature()
        nv = np.zeros(g.n_nodes)
        cnt = np.zeros(g.n_nodes)
        for e, (a, b) in enumerate(g.edges):
            nv[[a, b]] += c[e]
            cnt[[a, b]] += 1
        nv /= np.maximum(cnt, 1)
        label = "curvature"
    else:
        nv = model.node_means(X, values)
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(s, nv[order], color=CATEGORICAL[0], linewidth=2, marker="o", markersize=4)
    ax.set_xlabel("arc length along the filament", color=MUTED)
    ax.set_ylabel(label, color=MUTED)
    ax.grid(color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(colors=MUTED, labelsize=7)
    _title(ax, f"{label} along the filament")
    fig.tight_layout()
    return fig


def plot_cross_sections(model, X, *, n_sections=6, half_width=None, thickness=None, figsize=None):
    """Points around a filament seen in planes orthogonal to it at equally spaced nodes
    (co-moving orthonormal planes)."""
    plt = _plt()
    model = getattr(model, "agtm", model)
    from .crawling import ManifoldGraph

    X = np.asarray(X, float)
    g = ManifoldGraph(model.means_, model.graph.edges, model.graph.dim)
    order, s = chain_order(g)
    V = g.nodes[order]
    step = np.median(np.linalg.norm(np.diff(V, axis=0), axis=1))
    half_width = half_width or 3 * step
    thickness = thickness or step
    picks = np.unique(np.linspace(1, len(order) - 2, min(n_sections, max(len(order) - 2, 1))).astype(int))
    fig, axs = plt.subplots(1, len(picks), figsize=figsize or (2.8 * len(picks), 3.0), squeeze=False)
    prev = None
    for ax, i in zip(axs[0], picks):
        t = V[min(i + 1, len(V) - 1)] - V[max(i - 1, 0)]
        t /= np.linalg.norm(t)
        # orthonormal basis of the normal space, kept consistent between sections
        B = np.linalg.svd(np.eye(len(t)) - np.outer(t, t))[0][:, :len(t) - 1]
        if prev is not None:
            U, _, Wt = np.linalg.svd(B.T @ prev)
            B = B @ (U @ Wt)
        prev = B
        rel = X - V[i]
        m = (np.abs(rel @ t) < thickness / 2) & (np.linalg.norm(rel, axis=1) < half_width * 1.5)
        Q = rel[m] @ B[:, :2] if B.shape[1] >= 2 else np.c_[rel[m] @ B[:, 0], np.zeros(m.sum())]
        ax.scatter(*Q.T, s=3, color=CATEGORICAL[0], linewidths=0)
        ax.plot(0, 0, marker="+", color=INK, markersize=8)
        ax.set_xlim(-half_width, half_width)
        ax.set_ylim(-half_width, half_width)
        ax.set_aspect("equal")
        ax.grid(color=GRID, linewidth=0.5)
        ax.tick_params(colors=MUTED, labelsize=6)
        for sp in ax.spines.values():
            sp.set_color(GRID)
        _title(ax, f"s = {s[i]:.3g}")
    fig.suptitle("Cross-sections orthogonal to the filament (+ = node)", fontsize=9, color=INK, x=0.01, ha="left")
    fig.tight_layout()
    return fig


# ------------------------------------------------------------------ summary
def plot_summary(learner, *, axes=None, projection=None, names=None, max_points=40000, figsize=(15, 9)):
    """One figure with the outcome of every stage of a fitted :class:`MultiManifoldLearner`."""
    plt = _plt()
    L = learner
    X = L.X_
    v = View(X, axes, projection, names)
    fig = plt.figure(figsize=figsize)
    s = _sub(len(X), max_points)
    P = v(X[s])
    ax = v.new_axes(fig, 2, 3, 1)
    _scatter(ax, P, 0.4, color=CATEGORICAL[0], alpha=0.5)
    v.equal(ax, P)
    _title(ax, f"1. Data: {len(X):,} points")
    pre = L.preprocess_
    F = v(pre.diffused)
    ax = v.new_axes(fig, 2, 3, 2)
    sp = _sub(len(pre), max_points)
    _scatter(ax, v(pre.diffused[sp]), 0.5, color=CATEGORICAL[1])
    v.equal(ax, F, pad=1.05, robust=0.002)
    _title(ax, f"2. Preprocessed: {len(pre):,} points kept, diffused")
    ax = v.new_axes(fig, 2, 3, 3)
    idx = L.dim_result_.index[sp]
    for j in np.unique(idx):
        _scatter(ax, v(pre.diffused[sp][idx == j]), 0.6, color=DIM_COLORS.get(int(j), INK))
    v.equal(ax, F, pad=1.05, robust=0.002)
    _title(ax, "3. Dimensionality index")
    _legend(ax, [("undetermined" if j == 0 else f"{j}-D", DIM_COLORS.get(int(j), INK))
                 for j in np.unique(L.dim_result_.index)], loc="upper right")
    ax = v.new_axes(fig, 2, 3, 4)
    _scatter(ax, P, 0.3, color=NOISE, alpha=0.3)
    for k, m in enumerate(L.manifolds_):
        _edges(ax, v, m.graph, color_of(k), 0.8)
    v.equal(ax, F, pad=1.05, robust=0.002)
    _title(ax, f"4. Graphs: {len(L.manifolds_)} manifolds")
    ax = v.new_axes(fig, 2, 3, 5)
    lab = L.labels_[s]
    _scatter(ax, P[lab < 0], 0.3, color=NOISE, alpha=0.3)
    for k in range(len(L.manifolds_)):
        _scatter(ax, P[lab == k], 0.7, color=color_of(k))
    v.equal(ax, F, pad=1.05, robust=0.002)
    _title(ax, "5. Points attached to each manifold")
    ax = fig.add_subplot(2, 3, 6)
    for k, m in enumerate(L.manifolds_):
        if m.agtm is not None:
            ax.plot(m.agtm.loglik_, color=color_of(k), linewidth=1.4)
    ax.set_xlabel("EM iteration", color=MUTED)
    ax.set_ylabel("mean log-likelihood per point", color=MUTED)
    ax.grid(color=GRID, linewidth=0.5)
    for spn in ("top", "right"):
        ax.spines[spn].set_visible(False)
    ax.tick_params(colors=MUTED, labelsize=7)
    _title(ax, "6. AGTM training")
    fig.tight_layout()
    return fig
