"""Dimensionality index of points in a (diffused) point cloud (paper section 4.1)."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._util import as_points, build_tree, check_radius

__all__ = ["DimIndexResult", "dimensionality_index", "barycentric", "fisher_distance"]


@dataclass
class DimIndexResult:
    """Outcome of :func:`dimensionality_index`.

    ``index`` is the spatially smoothed soft index (the one used downstream), ``index_hard``
    the nearest-ideal-spectrum index. Both take values in 1..d, or 0 where the neighbourhood
    was too small to decide.
    """

    index: np.ndarray
    index_hard: np.ndarray
    proba: np.ndarray           # (n, d) soft index distribution P_i(j)
    proba_smooth: np.ndarray    # (n, d) spatially smoothed distribution P^S_i(j)
    eigenvalues: np.ndarray     # (n, d) normalised local spectrum, descending
    barycentric: np.ndarray     # (n, d) barycentric coordinates in the simplex S
    counts: np.ndarray          # (n,) neighbourhood sizes
    tangents: np.ndarray | None = None  # (n, n_vec, d) leading local eigenvectors

    def counts_per_dim(self) -> dict[int, int]:
        vals, cnt = np.unique(self.index, return_counts=True)
        return {int(v): int(c) for v, c in zip(vals, cnt)}


def fisher_distance(p, q) -> np.ndarray:
    """Fisher geodesic distance between points of the probability simplex (eq. 26)."""
    s = np.sum(np.sqrt(np.clip(p, 0, None) * np.clip(q, 0, None)), axis=-1)
    return 2.0 * np.arccos(np.clip(s, 0.0, 1.0))


def barycentric(spectrum) -> np.ndarray:
    """Barycentric coordinates of normalised, descending spectra with respect to the ideal
    spectra s_j = (1/j, ..., 1/j, 0, ..., 0):  alpha_j = j (lambda_j - lambda_{j+1})."""
    L = np.asarray(spectrum, dtype=float)
    nxt = np.concatenate([L[..., 1:], np.zeros(L.shape[:-1] + (1,))], axis=-1)
    return np.clip((L - nxt) * np.arange(1, L.shape[-1] + 1), 0.0, 1.0)


def _weighted_local_pca(tree, X, r, w, n_vec, chunk=2_000_000):
    """Local PCA in which neighbour ``j`` counts with weight ``w[j]`` -> counts, eigenvalues, vectors."""
    n, d = X.shape
    indptr, idx, _ = tree.radius_neighbors(X, float(r))
    counts = np.diff(indptr)
    ev = np.zeros((n, d))
    vec = np.zeros((n, n_vec, d))
    start = 0
    while start < n:                                    # rows in blocks of bounded size
        stop = max(int(np.searchsorted(indptr, indptr[start] + chunk, side="right")) - 1, start + 1)
        stop = min(stop, n)
        lo, hi = indptr[start], indptr[stop]
        j = idx[lo:hi]
        off = indptr[start:stop] - lo
        wj = w[j]
        sw = np.add.reduceat(wj, off)
        mean = np.add.reduceat(wj[:, None] * X[j], off, axis=0) / sw[:, None]
        Dc = X[j] - np.repeat(mean, counts[start:stop], axis=0)
        C = np.add.reduceat(wj[:, None, None] * Dc[:, :, None] * Dc[:, None, :], off, axis=0) / sw[:, None, None]
        val, V = np.linalg.eigh(C)
        ev[start:stop] = np.maximum(val[:, ::-1], 0.0)
        if n_vec:
            vec[start:stop] = np.transpose(V[:, :, ::-1], (0, 2, 1))[:, :n_vec]
        start = stop
    return counts, ev, vec


def _weighted_smooth(tree, X, P, r, sigma, valid, w, chunk=2_000_000):
    """Gaussian average of the rows of ``P`` over neighbours, neighbour ``j`` weighted by ``w[j]``."""
    n = len(X)
    indptr, idx, dist = tree.radius_neighbors(X, float(r))
    out = np.zeros_like(P)
    start = 0
    while start < n:
        stop = max(int(np.searchsorted(indptr, indptr[start] + chunk, side="right")) - 1, start + 1)
        stop = min(stop, n)
        lo, hi = indptr[start], indptr[stop]
        j = idx[lo:hi]
        g = np.exp(-dist[lo:hi] ** 2 / (2 * sigma ** 2)) * w[j] * valid[j]
        off = indptr[start:stop] - lo
        den = np.add.reduceat(g, off)
        out[start:stop] = np.add.reduceat(g[:, None] * P[j], off, axis=0) / np.maximum(den, 1e-300)[:, None]
        start = stop
    out[~valid.astype(bool)] = 0.0
    return out


def dimensionality_index(X, r: float, *, smooth_factor: float = 1.0, min_points: int | None = None,
                         n_vec: int = 0, tree=None, density_weight: bool = False) -> DimIndexResult:
    """Estimate the local intrinsic dimensionality of every point.

    Parameters
    ----------
    X : (n, d) array
        Points, normally the diffused set.
    r : float
        Radius of the neighbourhoods used for the local PCA.
    smooth_factor : float
        The soft index is smoothed over a ball of radius ``smooth_factor * r`` with a Gaussian
        of the same scale. Use 0 to disable smoothing.
    min_points : int, optional
        Minimum neighbourhood size (including the point itself) for a valid estimate.
        Default ``max(d + 1, 5)``.
    n_vec : int
        Number of leading local eigenvectors to return (``tangents``).
    density_weight : bool
        Weight every neighbour by the inverse of its own local density (the number of points in
        its ball of radius ``r``), in the local PCA and in the smoothing. The index then describes
        the shape of the region the data occupy and not where most of the points are: without
        it, a neighbourhood that straddles a dense and a sparse structure is described by the
        dense one alone, which moves the boundary between them into the sparse structure.
    """
    X = as_points(X)
    r = check_radius(r)
    n, d = X.shape
    if min_points is None:
        min_points = max(d + 1, 5)
    tree = tree if tree is not None else build_tree(X)
    if density_weight:
        w = 1.0 / np.maximum(tree.radius_count(X, r), 1)
        counts, ev, vec = _weighted_local_pca(tree, X, r, w, int(n_vec))
    else:
        counts, ev, vec = tree.local_pca(X, r, int(n_vec), int(min_points))
    tot = ev.sum(axis=1)
    valid = (counts >= min_points) & (tot > 0)
    L = np.zeros_like(ev)
    L[valid] = ev[valid] / tot[valid, None]

    # hard index: closest ideal spectrum s_j in the Fisher metric (eq. 27)
    js = np.arange(1, d + 1)
    csum = np.cumsum(np.sqrt(L), axis=1) / np.sqrt(js)
    dist_ideal = 2.0 * np.arccos(np.clip(csum, 0.0, 1.0))
    hard = np.where(valid, np.argmin(dist_ideal, axis=1) + 1, 0)

    # soft index: isotropic kernel on the vertices of the full simplex
    alpha = barycentric(L)
    kappa = 2.0 * np.arccos(1.0 / np.sqrt(d)) if d > 1 else 1.0
    dist_vertex = 2.0 * np.arccos(np.clip(np.sqrt(alpha), 0.0, 1.0))
    P = np.exp(-dist_vertex**2 / (2.0 * kappa**2))
    P /= P.sum(axis=1, keepdims=True)
    P[~valid] = 0.0

    if smooth_factor > 0:
        rs = smooth_factor * r
        PS = (_weighted_smooth(tree, X, P, rs, rs, valid.astype(float), w) if density_weight
              else tree.kernel_smooth(P, rs, rs, valid.astype(np.uint8)))
    else:
        PS = P.copy()
    index = np.where(valid, np.argmax(PS, axis=1) + 1, 0)
    return DimIndexResult(index=index, index_hard=hard, proba=P, proba_smooth=PS, eigenvalues=L,
                          barycentric=alpha, counts=counts, tangents=vec if n_vec > 0 else None)
