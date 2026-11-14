"""Two views of quantities carried by the points of a one-dimensional manifold (a filament or
stream), after Canducci et al., *1-DREAM: 1D Recovery, Extraction and Analysis of Manifolds in
noisy environments* (2022), section 4.

* **Bi-dimensional profile** (sec. 4.1): the mean of a quantity in cylindrical shells around
  every segment of the fitted curve, shown as one image with the position along the manifold on
  one axis and the distance from its core on the other.
* **Co-moving orthonormal frames** (sec. 4.2): a plane perpendicular to the curve, carrying a
  consistent pair of axes, is moved along a finely sampled spine; on every plane the quantities
  and the density of the fitted model are mapped on a square grid. The frames form a movie.

Both use the fitted 1-D :class:`manydream.AGTM` (the stream GTM): its trained mapping gives the spine
and its tangents. Means are Gaussian-weighted (sec. 4.3); the SPH weights of the paper are not
implemented. The functions work in any ambient dimension, with the axes of the paper in 3-D.
"""
from __future__ import annotations

import numpy as np

from ._util import as_points, build_tree

__all__ = ["spine", "bidimensional_profile", "comoving_frames", "plot_bidimensional_profile",
           "plot_comoving_frame", "animate_comoving_frames"]


def _model(model):
    model = getattr(model, "agtm", model)
    if model.graph.dim != 1:
        raise ValueError("these views are defined for one-dimensional manifolds")
    return model


def spine(model, upsample: int = 0) -> dict:
    """The fitted curve of a 1-D model, optionally refined between its centres.

    ``upsample`` latent points are inserted, equally spaced, in every interval between two
    consecutive centres and mapped to the data space with the trained mapping (eq. 53). Returns
    ``points`` (P, d), unit ``tangents`` (P, d), ``arc`` length (P,), ``is_node`` (P,) and
    ``order``, the node indices along the curve.
    """
    from .viz import chain_order

    model = _model(model)
    order, _ = chain_order(model.graph)
    edges = {tuple(e): k for k, e in enumerate(np.asarray(model.graph.edges))}
    pts, is_node = [], []
    for a, b in zip(order[:-1], order[1:]):
        s = np.arange(upsample + 1) / (upsample + 1.0)              # the first centre and the inserted points
        if (a, b) in edges:
            Phi = model.basis(edges[(a, b)], s)
        else:
            Phi = model.basis(edges[(b, a)], 1.0 - s)
        pts.append((model.W_ @ Phi).T)
        is_node += [True] + [False] * upsample
    pts.append(model.means_[order[-1:]])
    is_node.append(True)
    P = np.vstack(pts)
    seg = np.linalg.norm(np.diff(P, axis=0), axis=1)
    arc = np.r_[0.0, np.cumsum(seg)]
    T = np.gradient(P, arc, axis=0) if len(P) > 2 else np.repeat(np.diff(P, axis=0), len(P), axis=0)
    T /= np.maximum(np.linalg.norm(T, axis=1, keepdims=True), 1e-300)
    return {"points": P, "tangents": T, "arc": arc, "is_node": np.array(is_node), "order": order}


def _kernel_means(tree, values, Q, h, weights=None):
    """Gaussian-weighted mean of every column of ``values`` at the query points, and the weighted
    number density there. Query points with no data within 3 h get NaN (density 0)."""
    indptr, idx, dist = tree.radius_neighbors(np.ascontiguousarray(Q), 3.0 * h)
    w = np.exp(-dist ** 2 / (2 * h * h))
    if weights is not None:
        w = w * weights[idx]
    cnt = np.diff(indptr)
    row = np.repeat(np.arange(len(Q)), cnt)
    den = np.bincount(row, w, len(Q))
    out = np.full((len(Q), values.shape[1]), np.nan)
    ok = den > 0
    for k in range(values.shape[1]):
        v = values[idx, k]
        good = np.isfinite(v)
        num = np.bincount(row[good], w[good] * v[good], len(Q))
        dk = np.bincount(row[good], w[good], len(Q))
        out[:, k] = np.where(dk > 0, num / np.maximum(dk, 1e-300), np.nan)
    d = Q.shape[1]
    return out, np.where(ok, den / ((2 * np.pi) ** (d / 2) * h ** d), 0.0)


def _normal_basis(t, prev=None):
    """Orthonormal basis (d, d-1) of the space perpendicular to the unit vector ``t``. In 3-D it
    is the pair of eq. 55-56, which varies smoothly along a curve; otherwise (and where that pair
    degenerates) a basis rotated to stay as close as possible to the previous one."""
    d = len(t)
    if d == 3 and (t[0] ** 2 + t[1] ** 2) > 1e-8:
        u1 = np.array([t[1], -t[0], 0.0])
        u2 = np.array([t[0] * t[2], t[1] * t[2], -(t[0] ** 2 + t[1] ** 2)])
        return np.c_[u1 / np.linalg.norm(u1), u2 / np.linalg.norm(u2)]
    B = np.linalg.svd(np.eye(d) - np.outer(t, t))[0][:, :d - 1]
    if prev is not None and prev.shape == B.shape:
        U, _, Wt = np.linalg.svd(B.T @ prev)
        B = B @ (U @ Wt)
    return B


def _as_values(values, n):
    if values is None:
        return np.zeros((n, 0)), []
    if isinstance(values, dict):
        names = list(values)
        V = np.column_stack([np.asarray(values[k], float) for k in names])
    else:
        V = np.asarray(values, float)
        V = V[:, None] if V.ndim == 1 else V
        names = [f"value {k + 1}" for k in range(V.shape[1])]
    if len(V) != n:
        raise ValueError("values must hold one row per point")
    return V, names


def bidimensional_profile(model, X, values=None, *, r_max: float, n_radial: int = 19, n_samples: int = 600,
                          smoothing: float | None = None, weights=None, random_state=0) -> dict:
    """Bi-dimensional profile of a 1-D manifold (sec. 4.1).

    For every segment between two consecutive centres a cylinder of radius ``r_max`` around the
    segment, along the local tangent, is sampled uniformly (eq. 47-49). At every sample the
    Gaussian-weighted mean of each quantity over the data is computed, and the samples are
    averaged in ``n_radial - 1`` concentric shells (eq. 48, 50).

    Parameters
    ----------
    model : fitted 1-D AGTM (or a manifold of a learner)
    X : (n, d) points carrying the quantities
    values : (n,), (n, q) array or dict name -> (n,), optional
        Quantities to profile. The number density of the points is always returned.
    r_max : float
        Largest distance from the curve (the paper's r_M).
    n_radial : int
        Number of radial bin edges (the paper's c); there are ``n_radial - 1`` shells.
    n_samples : int
        Sample points per cylinder.
    smoothing : float, optional
        Scale of the Gaussian weights. Default: the width of a radial bin.
    weights : (n,) array, optional
        Weight of every point (e.g. its mass).

    Returns a dict with ``profiles`` (name -> (segments, shells) array, including ``"density"``),
    ``arc`` (segment edges along the curve) and ``radius`` (shell edges).
    """
    model = _model(model)
    X = as_points(X)
    V, names = _as_values(values, len(X))
    sp = spine(model, 0)
    P, T, arc = sp["points"], sp["tangents"], sp["arc"]
    d = X.shape[1]
    redges = np.linspace(0.0, float(r_max), int(n_radial))
    h = float(smoothing) if smoothing is not None else redges[1] - redges[0]
    rng = np.random.default_rng(random_state)
    # one uniform sample of the unit cylinder, reused for every segment (the paper's cloud C)
    u = rng.normal(size=(n_samples, d - 1))
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    disc = u * (rng.uniform(size=(n_samples, 1)) ** (1.0 / (d - 1))) * r_max
    axial = rng.uniform(size=n_samples)
    rad = np.linalg.norm(disc, axis=1)
    shell = np.clip(np.searchsorted(redges, rad, side="right") - 1, 0, n_radial - 2)
    tree = build_tree(X)
    w = None if weights is None else np.asarray(weights, float)
    nseg = len(P) - 1
    prof = {k: np.full((nseg, n_radial - 1), np.nan) for k in names + ["density"]}
    B = None
    for s in range(nseg):
        t = P[s + 1] - P[s]
        length = np.linalg.norm(t)
        t = t / max(length, 1e-300)
        B = _normal_basis(t, B)
        Q = P[s] + axial[:, None] * length * t + disc @ B.T          # eq. 49: scale, rotate, shift
        means, dens = _kernel_means(tree, V, Q, h, w)
        for i in range(n_radial - 1):
            m = shell == i
            if not m.any():
                continue
            prof["density"][s, i] = dens[m].mean()
            for k, nm in enumerate(names):
                col = means[m, k]
                if np.isfinite(col).any():
                    prof[nm][s, i] = np.nanmean(col)
    return {"profiles": prof, "arc": arc, "radius": redges, "tangents": T, "smoothing": h}


def comoving_frames(model, X, values=None, *, side: float, upsample: int = 5, n_pixels: int = 25,
                    depth: float | None = None, half_thickness: float | None = None,
                    smoothing: float | None = None, weights=None) -> dict:
    """Co-moving orthonormal coordinate frames along a 1-D manifold (sec. 4.2).

    The spine is refined with ``upsample`` points per interval (eq. 53). At every spine point a
    plane perpendicular to the tangent carries an ``n_pixels`` x ``n_pixels`` grid of side
    ``side`` centred on the spine (eq. 54-57). The points within ``half_thickness`` of the plane
    along the tangent and within ``depth`` of the spine across it are selected, and on every
    pixel the Gaussian-weighted mean of each quantity over them is computed, together with the
    density of the fitted model (eq. 58).

    Parameters
    ----------
    side : float
        Side of the plane (the paper's a).
    upsample : int
        Latent points inserted between consecutive centres (the paper's N_l = 5).
    n_pixels : int
        Pixels per side (the paper's M = 25); made odd so that the spine is on a pixel.
    depth : float, optional
        Largest distance from the spine of the points used (the paper's b). Default: half the
        diagonal of the plane.
    half_thickness : float, optional
        Half thickness of the slab of points around the plane. Default: half the spacing of the
        refined spine, as in the paper; increase it when the data are sparse.
    smoothing : float, optional
        Scale of the Gaussian weights. Default: 1.5 pixels.

    Returns a dict with, per frame, ``centers``, ``tangents``, ``axes`` (the two in-plane unit
    vectors), ``arc``, ``n_points``, the model density ``pdf`` (frames, M, M) and ``maps``
    (name -> (frames, M, M), including ``"density"``), plus the grid ``extent``.
    """
    model = _model(model)
    X = as_points(X)
    V, names = _as_values(values, len(X))
    M = int(n_pixels) | 1
    side = float(side)
    b = float(depth) if depth is not None else side * np.sqrt(2) / 2
    h = float(smoothing) if smoothing is not None else 1.5 * side / M
    sp = spine(model, int(upsample))
    P, T, arc = sp["points"], sp["tangents"], sp["arc"]
    half = float(half_thickness) if half_thickness is not None else 0.5 * float(np.median(np.diff(arc)))
    lin = (np.arange(M) - (M - 1) / 2) * (side / M)
    gu, gv = np.meshgrid(lin, lin, indexing="xy")
    tree = build_tree(X)
    w = None if weights is None else np.asarray(weights, float)
    reach = float(np.hypot(half, b))
    F = len(P)
    maps = {k: np.full((F, M, M), np.nan) for k in names + ["density"]}
    pdf = np.zeros((F, M, M))
    axes = np.zeros((F, 2, X.shape[1]))
    npts = np.zeros(F, dtype=int)
    B = None
    for f in range(F):
        B = _normal_basis(T[f], B)
        u1, u2 = B[:, 0], (B[:, 1] if B.shape[1] > 1 else np.zeros(len(T[f])))
        axes[f] = u1, u2
        grid = P[f] + gu.ravel()[:, None] * u1 + gv.ravel()[:, None] * u2
        pdf[f] = np.exp(model.score_samples(grid)).reshape(M, M)
        _, cand, _ = tree.radius_neighbors(P[f:f + 1], reach)
        rel = X[cand] - P[f]
        along = rel @ T[f]
        perp = rel - along[:, None] * T[f]
        sel = cand[(np.abs(along) <= half) & (np.linalg.norm(perp, axis=1) <= b)]
        npts[f] = len(sel)
        if len(sel) == 0:
            maps["density"][f] = 0.0
            continue
        sub = build_tree(np.ascontiguousarray(X[sel]))
        means, dens = _kernel_means(sub, V[sel], grid, h, None if w is None else w[sel])
        maps["density"][f] = dens.reshape(M, M)
        for k, nm in enumerate(names):
            maps[nm][f] = means[:, k].reshape(M, M)
    return {"centers": P, "tangents": T, "axes": axes, "arc": arc, "n_points": npts, "pdf": pdf, "maps": maps,
            "extent": (-side / 2, side / 2, -side / 2, side / 2), "side": side, "depth": b, "half_thickness": half,
            "smoothing": h, "is_node": sp["is_node"]}


# ------------------------------------------------------------------ figures
def plot_bidimensional_profile(result, *, names=None, log_density: bool = True, cmap="viridis", figsize=None):
    """Bi-dimensional profiles as images: position along the manifold against distance from its core."""
    from .viz import INK, MUTED, _plt

    plt = _plt()
    prof = result["profiles"]
    names = list(names) if names is not None else list(prof)
    fig, axs = plt.subplots(len(names), 1, figsize=figsize or (8.5, 1.9 * len(names) + 0.8), squeeze=False, sharex=True)
    for ax, nm in zip(axs[:, 0], names):
        Z = prof[nm].T
        label = nm
        if nm == "density" and log_density:
            Z = np.log10(np.where(Z > 0, Z, np.nan))
            label = "log10 density"
        im = ax.pcolormesh(result["arc"], result["radius"], Z, cmap=cmap, shading="flat")
        cb = fig.colorbar(im, ax=ax, pad=0.01)
        cb.set_label(label, fontsize=8, color=INK)
        cb.ax.tick_params(labelsize=7, colors=MUTED)
        ax.set_ylabel("distance from the core", fontsize=8, color=INK)
        ax.tick_params(labelsize=7, colors=MUTED)
    axs[-1, 0].set_xlabel("position along the manifold (arc length)", fontsize=8, color=INK)
    fig.tight_layout()
    return fig


def _limits(A, log=False):
    A = np.asarray(A, float)
    if log:
        A = np.log10(np.where(A > 0, A, np.nan))
    f = A[np.isfinite(A)]
    if not len(f):
        return 0.0, 1.0
    lo, hi = np.quantile(f, [0.02, 0.98])
    return (lo, hi) if hi > lo else (lo - 0.5, lo + 0.5)


def plot_comoving_frame(result, k: int, X=None, *, names=None, limits=None, cmap="viridis", fig=None, max_points=20000):
    """One co-moving frame: the quantities and the model density on the perpendicular plane, and
    where the plane is on the manifold and in the whole data set (needs ``X``)."""
    from .viz import INK, MUTED, NOISE, _plt

    plt = _plt()
    names = list(names) if names is not None else [n for n in result["maps"] if n != "density"][:2] + ["density"]
    panels = names + ["model density"]
    ncol = max(len(panels), 3 if X is not None else len(panels))
    fig = fig or plt.figure(figsize=(3.4 * ncol, 6.6))
    fig.clear()
    limits = limits or {}
    for j, nm in enumerate(panels):
        ax = fig.add_subplot(2, ncol, j + 1) if nm != "model density" else fig.add_subplot(2, ncol, ncol + 1)
        log = nm in ("density", "model density")
        Z = result["pdf"][k] if nm == "model density" else result["maps"][nm][k]
        Zp = np.log10(np.where(Z > 0, Z, np.nan)) if log else Z
        lo, hi = limits.get(nm) or _limits(result["pdf"] if nm == "model density" else result["maps"][nm], log)
        im = ax.imshow(Zp, origin="lower", extent=result["extent"], cmap=cmap, vmin=lo, vmax=hi)
        ax.plot(0, 0, marker="+", color="white", markersize=9, markeredgewidth=1.4)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.ax.tick_params(labelsize=7, colors=MUTED)
        ax.set_title(("log10 " if log else "") + nm, fontsize=9, color=INK)
        ax.set_xlabel("u1", fontsize=8, color=MUTED)
        ax.set_ylabel("u2", fontsize=8, color=MUTED)
        ax.tick_params(labelsize=7, colors=MUTED)
    if X is not None:
        X = np.asarray(X, float)
        c, P = result["centers"][k], result["centers"]
        d3 = min(X.shape[1], 3)
        pad = lambda A: np.c_[A[:, :d3], np.zeros((len(A), 3 - d3))]                 # noqa: E731
        near = np.linalg.norm(X[:, None, :] - P[None, ::max(len(P) // 60, 1), :], axis=2).min(1) < 1.5 * result["depth"] \
            if len(X) <= 40000 else np.linalg.norm(X - c, axis=1) < 6 * result["depth"]
        sub = np.random.default_rng(0).choice(len(X), min(len(X), max_points), replace=False)
        for col, pts, title in ((ncol + 2, X[near], "position on the manifold"), (ncol + 3, X[sub], "position in the data set")):
            if col > 2 * ncol:
                break
            ax = fig.add_subplot(2, ncol, col, projection="3d", computed_zorder=False)
            Q = pad(pts)
            ax.scatter(*Q.T, s=1.5, color=NOISE, linewidths=0, depthshade=False, alpha=0.45, zorder=1)
            ax.plot(*pad(P).T, color=INK, linewidth=1.2, zorder=2)
            ax.scatter(*pad(c[None]).T, s=90, color="#e34948", edgecolors="white", linewidths=0.8, depthshade=False, zorder=4)
            if title == "position on the manifold":
                a = result["side"] / 2
                u1, u2 = result["axes"][k]
                sq = np.array([c + s1 * a * u1 + s2 * a * u2 for s1, s2 in ((-1, -1), (1, -1), (1, 1), (-1, 1), (-1, -1))])
                ax.plot(*pad(sq).T, color="#2a78d6", linewidth=1.4, zorder=3)
            ax.set_title(title, fontsize=9, color=INK)
            ax.set_axis_off()
            ax.set_box_aspect((1, 1, 1))
    fig.suptitle(f"frame {k + 1} of {len(result['centers'])} · arc length {result['arc'][k]:.3g} · "
                 f"{result['n_points'][k]} points on the plane", fontsize=10, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig


def animate_comoving_frames(result, path, X=None, *, names=None, fps: int = 10, cmap="viridis", dpi: int = 100):
    """Join the co-moving frames into a movie (``.mp4`` with ffmpeg, otherwise ``.gif``).
    The colour scales are fixed over the whole movie. Returns the path written."""
    from pathlib import Path

    from matplotlib import animation

    from .viz import _plt

    plt = _plt()
    path = Path(path)
    names_ = list(names) if names is not None else [n for n in result["maps"] if n != "density"][:2] + ["density"]
    limits = {nm: _limits(result["maps"][nm], nm == "density") for nm in names_}
    limits["model density"] = _limits(result["pdf"], True)
    use_ffmpeg = path.suffix.lower() == ".mp4" and animation.writers.is_available("ffmpeg")
    if path.suffix.lower() == ".mp4" and not use_ffmpeg:
        path = path.with_suffix(".gif")
    writer = (animation.FFMpegWriter(fps=fps, codec="libx264", extra_args=["-pix_fmt", "yuv420p"]) if use_ffmpeg
              else animation.PillowWriter(fps=fps))
    ncol = max(len(names_) + 1, 3 if X is not None else len(names_) + 1)
    fig = plt.figure(figsize=(3.4 * ncol, 6.6))
    with writer.saving(fig, str(path), dpi=dpi):
        for k in range(len(result["centers"])):
            plot_comoving_frame(result, k, X, names=names_, limits=limits, cmap=cmap, fig=fig)
            writer.grab_frame()
    plt.close(fig)
    return path
