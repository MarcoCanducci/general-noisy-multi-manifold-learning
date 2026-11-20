"""Exact LAAT and its Markov chain, after Taghribi et al. (IEEE TKDE 35(6), 2023), equations (1)-(12).

This module follows the paper literally, with a fixed neighbourhood radius: alignment weights
(2), normalised eigenvalues (3), preference (4), relative preference (5), relative pheromone
(6), movement preference (7), Boltzmann jump probabilities (8), deposit (10) and evaporation
(1). It is meant for studying the method; for preprocessing data use :class:`manydream.LAAT`,
which is faster and differs from the paper in a few documented places.

Without reinforcement (``kappa = 1``) the ants form a Markov chain (11) whose stationary
distribution (12) is available from :meth:`LAATGraph.stationary` without simulating any ant.

Two choices that the paper leaves fixed can be varied here:

* the normalisation within a neighbourhood: ``norm="sum"`` (the paper) or ``norm="max"``.
  Dividing by the maximum makes the inverse temperature ``beta`` act alike in sparse and dense
  neighbourhoods; with the sum, ``beta`` is effectively divided by the number of neighbours;
* where the ants start (:meth:`LAATGraph.start_distribution`) and how long each walks
  (``length`` in :func:`run`).

Points are rows.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.sparse import csgraph
from scipy.sparse.linalg import spsolve
from scipy.spatial import cKDTree

__all__ = ["LAATGraph", "run", "walk_lengths", "START_RULES"]

START_RULES = ("uniform", "median", "density", "stationary", "complement", "inverse", "one_minus_pi")


class LAATGraph:
    """Neighbourhoods of a data set at radius ``r`` and the quantities that do not depend on the pheromone.

    Parameters
    ----------
    X : (n, D) array
    r : float
        Neighbourhood radius.
    min_nb : int, optional
        Points with fewer neighbours are removed from the graph (repeatedly, until every
        remaining point has enough remaining neighbours). Default ``D + 1``.
    norm : {"sum", "max"}
        Normalisation of the preference and of the pheromone within a neighbourhood.

    Attributes
    ----------
    valid : (n,) bool
        Points that take part in the walks.
    row, col : (E,) int
        The directed edges ``row -> col`` between valid neighbours, sorted by ``row``.
    deg : (n,) int
        Number of neighbours of every point.
    E, Ebar : (E,) float
        Preference (4) and relative preference (5) of every edge.
    """

    def __init__(self, X, r: float, min_nb: int | None = None, norm: str = "sum"):
        if norm not in ("sum", "max"):
            raise ValueError("norm must be 'sum' or 'max'")
        self.norm = norm
        X = np.asarray(X, float)
        n, D = X.shape
        min_nb = D + 1 if min_nb is None else min_nb
        pairs = cKDTree(X).query_pairs(r, output_type="ndarray")
        row = np.r_[pairs[:, 0], pairs[:, 1]]
        col = np.r_[pairs[:, 1], pairs[:, 0]]
        valid = np.ones(n, bool)
        while True:                                   # drop points with too few valid neighbours, until stable
            keep = valid[row] & valid[col]
            deg = np.bincount(row[keep], minlength=n)
            new = deg >= min_nb
            if (new == valid).all():
                break
            valid = new
        row, col = row[keep], col[keep]
        order = np.lexsort((col, row))
        row, col = row[order], col[order]
        deg = np.bincount(row, minlength=n)
        indptr = np.r_[0, np.cumsum(deg)]
        self.X, self.n, self.D, self.r = X, n, D, r
        self.valid, self.row, self.col, self.deg, self.indptr = valid, row, col, deg, indptr
        self.first = indptr[:-1]
        self.last = indptr[1:] - 1
        self._starts = self.first[deg > 0]
        self._counts = deg[deg > 0]

        # local PCA of every valid neighbourhood (the point and its neighbours), equations (2)-(4)
        E = np.zeros(len(row))
        for i in np.flatnonzero(valid):
            a, b = indptr[i], indptr[i + 1]
            lam, V = np.linalg.eigh(np.cov(np.vstack([X[col[a:b]], X[i]]).T))
            lam = np.maximum(lam, 0)
            d = X[col[a:b]] - X[i]
            cos = np.abs(d @ V) / np.linalg.norm(d, axis=1, keepdims=True)
            E[a:b] = (cos / cos.sum(1, keepdims=True)) @ (lam / lam.sum())
        self.E = E
        self.Ebar = E / self._reduce(E)                # equation (5)

    def _reduce(self, v):
        """Sum (or maximum) of ``v`` over the edges of each point, repeated on those edges."""
        fn = np.add.reduceat if self.norm == "sum" else np.maximum.reduceat
        return fn(v, self._starts).repeat(self._counts)

    def _row_sum(self, v):
        return np.add.reduceat(v, self._starts).repeat(self._counts)

    def probabilities(self, F, kappa: float, beta: float) -> np.ndarray:
        """Jump probabilities (8) on every edge, given the pheromone ``F`` on the points."""
        Fj = F[self.col]
        Fbar = Fj / self._reduce(Fj)                                   # equation (6)
        V = (1 - kappa) * Fbar + kappa * self.Ebar                     # equation (7)
        vmax = np.maximum.reduceat(V, self._starts).repeat(self._counts)
        p = np.exp(beta * (V - vmax))
        return p / self._row_sum(p)

    def transition_matrix(self, beta: float) -> sparse.csr_matrix:
        """Transition matrix of the Markov chain without reinforcement (``kappa = 1``), equation (11)."""
        p = self.probabilities(np.ones(self.n), 1.0, beta)
        return sparse.csr_matrix((p, (self.row, self.col)), shape=(self.n, self.n))

    def stationary(self, beta: float, weights: str = "size") -> np.ndarray:
        """Stationary distribution of the ``kappa = 1`` chain, equation (12), by a sparse linear solve.

        The chain is ergodic only within a connected component of the neighbourhood graph.
        The distributions of the components are combined with weights proportional to their
        number of points (``weights="size"``) or equal (``"equal"``). Points outside the graph
        get zero. The number of components is stored in ``n_components``.
        """
        P = self.transition_matrix(beta)
        ncomp, lab = csgraph.connected_components(P, directed=False)
        comps = [np.flatnonzero(lab == c) for c in range(ncomp)]
        comps = [c for c in comps if self.valid[c].all() and len(c) > 1]
        pi = np.zeros(self.n)
        for c in comps:
            m = len(c)
            A = (P[c][:, c].tocsc().T - sparse.identity(m, format="csc")).tolil()
            A[m - 1, :] = 1.0                           # replace one equation by sum(pi) = 1
            b = np.zeros(m)
            b[m - 1] = 1.0
            pc = np.maximum(spsolve(A.tocsc(), b), 0)
            pi[c] = pc / pc.sum() * (m if weights == "size" else 1.0)
        self.n_components = len(comps)
        self.P1 = P
        return pi / pi.sum()

    def start_distribution(self, rule: str, beta: float | None = None) -> np.ndarray:
        """Probability of every point to be the start of an ant.

        ``"uniform"``: every valid point; ``"median"``: points with at least the median number
        of neighbours (equation 9); ``"density"``: proportional to the number of neighbours;
        ``"stationary"``: the stationary distribution ``pi``; ``"complement"``: proportional
        to ``max(pi) - pi``; ``"inverse"``: proportional to ``1 / pi``; ``"one_minus_pi"``:
        proportional to ``1 - pi`` (nearly uniform, since every ``pi`` is about ``1/n``).
        The last four need ``beta``.
        """
        v = self.valid.astype(float)
        if rule == "uniform":
            w = v
        elif rule == "median":
            w = v * (self.deg >= np.median(self.deg))
        elif rule == "density":
            w = v * self.deg
        elif rule in ("stationary", "complement", "inverse", "one_minus_pi"):
            if beta is None:
                raise ValueError(f"start rule '{rule}' needs beta")
            pi = self.stationary(beta)
            if rule == "stationary":
                w = pi
            elif rule == "complement":
                w = v * (pi.max() - pi)
            elif rule == "inverse":
                w = np.where(pi > 0, 1.0 / np.maximum(pi, 1e-300), 0.0)
            else:
                w = v * (1.0 - pi)
        else:
            raise ValueError(f"unknown start rule '{rule}'; choose from {START_RULES}")
        return w / w.sum()


def walk_lengths(weight, total: int, lo: float = 0.1, hi: float = 5.0) -> np.ndarray:
    """Whole numbers of jumps proportional to ``weight``, clipped to ``[lo, hi]`` times the mean and summing to about ``total``."""
    w = np.asarray(weight, float)
    m = total / len(w)
    L = w / w.mean() * m
    for _ in range(20):                                 # clip, then rescale the unclipped ones to keep the total
        L = np.clip(L, lo * m, hi * m)
        free = (L > lo * m) & (L < hi * m)
        if not free.any() or abs(L.sum() - total) < 1:
            break
        L[free] *= max((total - L[~free].sum()) / L[free].sum(), 1e-9)
    return np.maximum(np.round(np.clip(L, lo * m, hi * m)).astype(int), 1)


def run(G: LAATGraph, start: str | None = "median", kappa: float = 0.8, beta: float = 100.0, n_ants: int = 100,
        n_steps: int = 200, n_epochs: int = 30, zeta: float = 0.1, phi: float | None = None, seed=0,
        start_p=None, record=None, length=None) -> np.ndarray:
    """Run LAAT and return the pheromone on every point after the last epoch.

    Parameters
    ----------
    G : LAATGraph
    start : str
        Start rule (see :meth:`LAATGraph.start_distribution`); ignored if ``start_p`` is given.
    kappa : float
        Weight of the alignment against the pheromone in (7); 1 gives the Markov chain.
    beta : float
        Inverse temperature of the jump probabilities (8).
    n_ants, n_steps, n_epochs : int
        Ants per epoch, jumps per ant and epoch, and epochs.
    zeta : float
        Evaporation rate (1).
    phi : float, optional
        Pheromone deposited per visit (10). Default ``n / (n_ants * n_steps)``: about one unit
        per point and epoch.
    seed : int
    start_p : (n,) array, optional
        Start probabilities, instead of a named rule.
    record : callable, optional
        Called as ``record(epoch, F)`` after every epoch.
    length : callable, optional
        ``length(start_points) -> weights``: gives every ant its own number of jumps,
        proportional to the weight (through :func:`walk_lengths`), with the same total
        ``n_ants * n_steps``.
    """
    rng = np.random.default_rng(seed)
    n = G.n
    F = np.ones(n)                                      # F(0) = 1
    phi = n / float(n_ants * n_steps) if phi is None else phi
    sp = G.start_distribution(start, beta) if start_p is None else np.asarray(start_p, float)
    cs = np.cumsum(sp)
    for ep in range(n_epochs):
        p = G.probabilities(F, kappa, beta)
        cum = np.cumsum(p)
        base = np.where(G.first > 0, cum[np.maximum(G.first - 1, 0)], 0.0)
        ants = np.minimum(np.searchsorted(cs, rng.random(n_ants), side="right"), n - 1)
        visits = np.zeros(n)
        if length is None:
            for _ in range(n_steps):
                e = np.searchsorted(cum, base[ants] + rng.random(n_ants) * (cum[G.last[ants]] - base[ants]), side="right")
                ants = G.col[np.clip(e, G.first[ants], G.last[ants])]
                np.add.at(visits, ants, 1.0)
        else:                                           # every ant has its own number of jumps
            L = walk_lengths(length(ants), n_ants * n_steps)
            order = np.argsort(-L)
            ants, L = ants[order], L[order]
            for t in range(int(L[0])):
                k = int(np.searchsorted(-L, -t, side="left"))        # the ants that still have jumps left
                a = ants[:k]
                e = np.searchsorted(cum, base[a] + rng.random(k) * (cum[G.last[a]] - base[a]), side="right")
                ants[:k] = G.col[np.clip(e, G.first[a], G.last[a])]
                np.add.at(visits, ants[:k], 1.0)
        F = (1 - zeta) * (F + phi * visits)             # deposit (10), then evaporation (1)
        if record is not None:
            record(ep, F)
    return F
