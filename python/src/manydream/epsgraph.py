"""Efficient alternative to manifold crawling (paper section 4.3)."""
from __future__ import annotations

import numpy as np

from ._util import as_points, build_tree, check_radius, rng_from
from .crawling import ManifoldGraph, attach_points

__all__ = ["epsilon_graph"]


def epsilon_graph(diffused, r: float, dim: int, *, noisy=None, angle_threshold: float = 2.0 / 3.0,
                  min_points: int = 5, min_nodes: int | None = None, random_state=None) -> list[ManifoldGraph]:
    """Build manifold graphs from a ball cover of the diffused points.

    The points are covered by balls of radius ``r / 2``; ball means become nodes. Two nodes
    are connected when they are at most ``r`` apart and the principal angle between their
    local tangent spaces, ``arccos |det(U_m^T U_n)|``, is below ``angle_threshold * pi / 2``.
    Connected components of the resulting graph are returned as separate manifolds.
    """
    D = as_points(diffused)
    r = check_radius(r)
    n, d = D.shape
    if not 1 <= dim <= d:
        raise ValueError(f"dim must be between 1 and {d}")
    if min_nodes is None:
        min_nodes = 4 * dim + 4
    tree = build_tree(D)
    labels, _ = tree.ball_cover(r / 2.0, rng_from(random_state).permutation(n).astype(np.int64))
    nb = int(labels.max()) + 1
    cnt = np.bincount(labels, minlength=nb)
    centres = np.stack([np.bincount(labels, weights=D[:, k], minlength=nb) for k in range(d)], 1) / cnt[:, None]
    # local tangent spaces: PCA of the points within r/2 of each centre
    counts, ev, vec = tree.local_pca(centres, r / 2.0, dim, max(min_points, dim + 1))
    ok = (counts >= max(min_points, dim + 1)) & (ev[:, dim - 1] > 0)
    centres, vec = np.ascontiguousarray(centres[ok]), vec[ok]
    if len(centres) == 0:
        return []
    ct = build_tree(centres)
    indptr, ind, _ = ct.radius_neighbors(centres, r)
    src = np.repeat(np.arange(len(centres)), np.diff(indptr))
    keep = src < ind
    a, b = src[keep], ind[keep]
    M = np.einsum("eik,ejk->eij", vec[a], vec[b])
    ang = np.arccos(np.clip(np.abs(np.linalg.det(M)), 0.0, 1.0))
    good = ang < angle_threshold * np.pi / 2.0
    full = ManifoldGraph(centres, np.column_stack([a[good], b[good]]), dim, vec)
    graphs = [g for g in full.components() if g.n_nodes >= min_nodes]
    graphs.sort(key=lambda g: -g.n_nodes)
    if noisy is not None:
        Q = as_points(noisy, "noisy")
        qt = build_tree(Q)
        for g in graphs:
            g.points = attach_points(g, Q, r, qt)
    return graphs
