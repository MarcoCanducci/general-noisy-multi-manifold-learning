"""Abstract GTM: manifold-aligned density model on a graph-structured latent space (paper section 3.1)."""
from __future__ import annotations

import warnings

import numpy as np
from scipy import linalg, sparse
from scipy.sparse import csgraph
from scipy.sparse.linalg import LinearOperator, cg
from scipy.special import logsumexp

from . import _core
from ._util import as_points, build_tree, rng_from
from .crawling import ManifoldGraph

__all__ = ["AGTM", "MultiManifoldModel", "DegenerateModelWarning", "DegenerateModelError", "epsilon_net"]


class DegenerateModelWarning(UserWarning):
    """A mixture component became degenerate during EM."""


class DegenerateModelError(RuntimeError):
    """Raised when ``on_degenerate='raise'`` and a component degenerates."""


def epsilon_net(adjacency, epsilon: int = 1) -> np.ndarray:
    """Greedy epsilon-net of a graph: centres are more than ``epsilon`` hops apart and every
    node is within ``epsilon`` hops of a centre. Nodes are visited in breadth-first order."""
    A = sparse.csr_matrix(adjacency)
    n = A.shape[0]
    covered = np.zeros(n, bool)
    centres = []
    order = []
    seen = np.zeros(n, bool)
    for s in range(n):
        if not seen[s]:
            o = csgraph.breadth_first_order(A, s, directed=False, return_predecessors=False)
            seen[o] = True
            order.extend(o.tolist())
    indptr, indices = A.indptr, A.indices
    stamp = np.full(n, -1, dtype=np.int64)
    for v in order:
        if covered[v]:
            continue
        centres.append(v)
        frontier = [v]
        covered[v] = True
        stamp[v] = v
        for _ in range(int(epsilon)):
            nxt = []
            for u in frontier:
                for w in indices[indptr[u]:indptr[u + 1]]:
                    if stamp[w] != v:
                        stamp[w] = v
                        covered[w] = True
                        nxt.append(w)
            frontier = nxt
    return np.asarray(centres, dtype=np.int64)


class AGTM:
    """Abstract Generative Topographic Mapping of one noisy manifold.

    The latent space is the abstract graph of a :class:`ManifoldGraph`. Its nodes are mapped
    to data space by ``f(v) = W phi(v)`` with radial basis functions of the graph distance,
    and each mapped node carries a Gaussian whose covariance ``zeta_i * Sigma_i`` is aligned
    with the manifold. ``W`` and ``zeta`` are fitted by EM.

    Parameters
    ----------
    graph : ManifoldGraph
    epsilon : int
        Spacing (in hops) of the epsilon-net of basis-function centres.
    gamma : float, optional
        Basis-function width in hops; default ``2 * epsilon``.
    normalize : bool
        Normalise the basis functions to sum to one at every node (default). On the irregular
        graphs produced by crawling, plain radial basis functions cannot reproduce even the
        node positions; the normalised ones can.
    ridge : float
        Regularisation of the weights (chi in eq. 17).
    max_iter, tol :
        EM stops after ``max_iter`` iterations or when the mean log-likelihood per point
        improves by less than ``tol``.
    cov_floor : float
        Eigenvalues of each local covariance are floored at ``cov_floor`` times its largest one.
    zeta_bounds : (lo, hi), optional
        Clamp for the covariance scales. Unbounded by default.
    on_degenerate : {"freeze", "stop", "raise"}
        What to do when a component degenerates (see ``diagnostics_``).
    min_resp : float, optional
        Minimum effective number of points per component; default ``d + 1``.
    collapse_ratio, explode_ratio : float
        A component whose scale ``zeta_i`` leaves ``[collapse_ratio, explode_ratio]`` (relative to
        its local data covariance) is treated as degenerate: it is collapsing onto a few points
        or drifting away from the data.
    weighted_geodesics : bool
        Use embedded edge lengths (in units of their mean) instead of hop counts.
    prune : float or None
        Components whose log density at a point is more than ``prune`` below their peak are
        skipped in the E-step (None evaluates every pair).
    solver : {"auto", "direct", "cg"}
        Linear solver for the W update.
    """

    def __init__(self, graph: ManifoldGraph, *, epsilon: int = 1, gamma: float | None = None,
                 ridge: float = 1e-5, max_iter: int = 100, tol: float = 5e-4, cov_floor: float = 1e-3,
                 cov_scale: float | str = "edge", zeta_bounds=None, on_degenerate: str = "freeze",
                 min_resp: float | None = None, weighted_geodesics: bool = False, solver: str = "auto",
                 basis_cutoff: float = 1e-6, normalize: bool = True, prune: float | None = 30.0,
                 collapse_ratio: float = 1e-2, explode_ratio: float = 1e2):
        if on_degenerate not in ("freeze", "stop", "raise"):
            raise ValueError("on_degenerate must be 'freeze', 'stop' or 'raise'")
        if graph.n_nodes < 2 or graph.n_edges < 1:
            raise ValueError("the graph needs at least two connected nodes")
        self.graph = graph
        self.epsilon = int(epsilon)
        self.gamma = float(gamma) if gamma is not None else 2.0 * self.epsilon
        self.ridge = ridge
        self.max_iter = max_iter
        self.tol = tol
        self.cov_floor = cov_floor
        self.cov_scale = cov_scale
        self.zeta_bounds = zeta_bounds
        self.on_degenerate = on_degenerate
        self.min_resp = min_resp
        self.weighted_geodesics = weighted_geodesics
        self.solver = solver
        self.basis_cutoff = basis_cutoff
        self.normalize = normalize
        self.prune = prune
        self.collapse_ratio = collapse_ratio
        self.explode_ratio = explode_ratio
        self._build_basis()

    # ------------------------------------------------------------------ latent space
    def _build_basis(self):
        g = self.graph
        A = g.adjacency(weighted=self.weighted_geodesics)
        if self.weighted_geodesics:
            A = A / g.edge_lengths.mean()
        self.centres_ = epsilon_net(g.adjacency(), self.epsilon)
        cut = self.gamma * np.sqrt(-np.log(self.basis_cutoff))
        D = csgraph.dijkstra(A, directed=False, indices=self.centres_, unweighted=not self.weighted_geodesics,
                             limit=cut)
        self._geo = D                                   # (M, K) graph distances centre -> node
        Phi = np.exp(-(D / self.gamma) ** 2)
        Phi[~np.isfinite(D)] = 0.0
        if self.normalize:
            Phi /= Phi.sum(0, keepdims=True)
        self.Phi_ = Phi                                  # (M, K)
        self._Phi_sp = sparse.csr_matrix(Phi)

    def basis(self, edge: int, s) -> np.ndarray:
        """Basis functions at fraction(s) ``s`` in [0, 1] along edge number ``edge`` -> (M, len(s))."""
        a, b = self.graph.edges[edge]
        s = np.atleast_1d(np.asarray(s, float))
        w = 1.0
        if self.weighted_geodesics:
            w = self.graph.edge_lengths[edge] / self.graph.edge_lengths.mean()
        D = np.minimum(self._geo[:, [a]] + w * s[None, :], self._geo[:, [b]] + w * (1.0 - s[None, :]))
        P = np.exp(-(D / self.gamma) ** 2)
        P[~np.isfinite(D)] = 0.0
        if self.normalize:
            P /= P.sum(0, keepdims=True)
        return P

    # ------------------------------------------------------------------ fitting
    def _init_weights(self):
        Phi, V = self.Phi_, self.graph.nodes
        M = Phi.shape[0]
        G = Phi @ Phi.T + self.ridge * np.eye(M)
        self.W_ = linalg.solve(G, Phi @ V, assume_a="pos").T      # (d, M)
        self.means_ = (self.W_ @ Phi).T                            # (K, d)

    def _local_covariances(self, T, tree, weights=None):
        g = self.graph
        V = self.means_
        K, d = V.shape
        if self.cov_scale == "edge":
            el = np.linalg.norm(V[g.edges[:, 0]] - V[g.edges[:, 1]], axis=1)
            deg = np.maximum(g.degree, 1)
            ell = (np.bincount(g.edges[:, 0], el, K) + np.bincount(g.edges[:, 1], el, K)) / deg
            ref = float(g.edge_lengths.mean())          # scale of the crawled graph itself
            ell[g.degree == 0] = ref
            ell = np.maximum(ell, 0.1 * ref)
        else:
            ell = np.full(K, float(self.cov_scale))
        S, mass = tree.weighted_cov(V, ell, 3.0, 5.0, weights)
        # regularise: isotropic fallback where no data, relative eigenvalue floor elsewhere
        bad = ~np.isfinite(S).all(axis=(1, 2))
        S[bad] = 0.0
        w, U = np.linalg.eigh(S)
        top = w[:, -1]
        empty = ~(top > 0) | (mass < 1e-12) | bad
        scale = np.where(empty, ell**2 / 3.0, top)
        w = np.where(empty[:, None], scale[:, None], np.maximum(w, self.cov_floor * scale[:, None]))
        self.Sigma_ = np.einsum("kab,kb,kcb->kac", U, w, U)
        self._Sigma_inv = np.einsum("kab,kb,kcb->kac", U, 1.0 / w, U)
        self._logdet_Sigma = np.log(w).sum(1)
        self._wmax = w.max(1)
        self.cov_scale_ = ell
        return empty

    def _gauss_params(self, zeta):
        K, d = self.means_.shape
        L = np.linalg.cholesky(self.Sigma_ * zeta[:, None, None])
        Wh = np.linalg.inv(L)
        logdet = self._logdet_Sigma + d * np.log(zeta)
        logc = -np.log(K) - 0.5 * logdet - 0.5 * d * np.log(2 * np.pi)
        radii = None
        if self.prune:
            # beyond this distance a component's log density is more than `prune` below its peak
            radii = np.sqrt(2.0 * self.prune * zeta * self._wmax)
        return Wh, logc, radii

    def _solve_W(self, G, S, zeta):
        d, M = self.W_.shape
        K = len(G)
        Lam = self._Sigma_inv / zeta[:, None, None]                  # (K, d, d)
        rhs = np.einsum("kab,kb->ka", Lam, S)                         # (K, d)
        B = (self._Phi_sp @ rhs).T                                    # (d, M)
        # purely numerical regularisation, relative to the typical diagonal entry of the system
        lam = 1e-9 * max(float((G * np.einsum("kaa->k", Lam)).mean()) / d, 1e-300)
        use_direct = self.solver == "direct" or (self.solver == "auto" and d * M <= 1000)
        if use_direct:
            Phi = self.Phi_
            A = np.empty((M, d, M, d))
            for a in range(d):
                for b in range(a, d):
                    blk = (Phi * (G * Lam[:, a, b])) @ Phi.T
                    A[:, a, :, b] = blk
                    A[:, b, :, a] = blk
            A = A.reshape(M * d, M * d)
            A[np.diag_indices_from(A)] += lam
            b = B.T.reshape(-1)
            try:
                x = linalg.solve(A, b, assume_a="pos")
            except linalg.LinAlgError:
                x = linalg.lstsq(A, b)[0]
            return x.reshape(M, d).T
        Psp = self._Phi_sp
        GL = Lam * G[:, None, None]

        def matvec(x):
            W = x.reshape(M, d)                                       # rows = basis functions
            Vk = Psp.T @ W                                            # (K, d)
            Y = np.einsum("kab,kb->ka", GL, Vk)
            return (Psp @ Y + lam * W).reshape(-1)

        op = LinearOperator((M * d, M * d), matvec=matvec, dtype=float)
        x, info = cg(op, B.T.reshape(-1), x0=self.W_.T.reshape(-1), rtol=1e-5, maxiter=300)
        return x.reshape(M, d).T

    # ---- stepwise training interface (used by fit and by HierarchicalAGTM)
    def begin(self, T, sample_weight=None) -> "AGTM":
        """Initialise the model on data ``T`` without running EM.

        ``sample_weight`` (n,) weights every point's contribution (used for hierarchical
        training and for undoing a thinning of the data).
        """
        T = as_points(T, "T")
        N, d = T.shape
        if d != self.graph.nodes.shape[1]:
            raise ValueError("data and graph have different ambient dimensions")
        w = None
        if sample_weight is not None:
            w = np.ascontiguousarray(sample_weight, dtype=float)
            if w.shape != (N,) or (w < 0).any() or not np.isfinite(w).all() or w.sum() <= 0:
                raise ValueError("sample_weight must hold one non-negative finite weight per point")
        K = self.graph.n_nodes
        self._T, self._w = T, w
        self._min_resp = self.min_resp if self.min_resp is not None else d + 1.0
        self._init_weights()
        empty = self._local_covariances(T, build_tree(T), w)
        self._zeta = np.ones(K)
        self._frozen = np.zeros(K, bool)
        self._reasons = {int(k): "no data near node at initialisation" for k in np.flatnonzero(empty)}
        self._events, self._history, self._state, self._status, self._G = [], [], None, "max_iter", None
        self._sync()
        return self

    def _sync(self):
        """Refresh the quantities derived from the current parameters."""
        self.zeta_ = self._zeta
        self.covariances_ = self.Sigma_ * self._zeta[:, None, None]
        self._Wh, self._logc, self._radii = self._gauss_params(self._zeta)

    def _flag(self, ks, why, it):
        new = [int(k) for k in ks if int(k) not in self._reasons or not self._frozen[k]]
        for k in ks:
            self._reasons.setdefault(int(k), why)
        if new:
            self._events.append({"iteration": it, "components": new, "reason": why})

    def _estep(self, sample_weight=None):
        """Per-point log density under the current parameters and the (weighted) statistics."""
        w = self._w if sample_weight is None else sample_weight
        logp, G, S, C, _ = _core.agtm_estep(self._T, self.means_, self._Wh, self._logc, None, True, self._radii, w)
        self._G = G
        return logp, G, S, C

    def _mstep(self, G, S, C, it: int = 0) -> bool:
        """Update weights (eq. 9) and scales (eq. 10). Returns False if training must stop."""
        d = self._T.shape[1]
        K = len(G)
        zeta = self._zeta
        Gs = np.maximum(G, 1e-300)
        W_new = self._solve_W(G, S, zeta)
        V = (W_new @ self.Phi_).T
        Mk = C - np.einsum("ka,kb->kab", S, V) - np.einsum("ka,kb->kab", V, S) \
            + G[:, None, None] * np.einsum("ka,kb->kab", V, V)
        z_new = np.einsum("kab,kba->k", self._Sigma_inv, Mk) / (d * Gs)
        starved = G < self._min_resp
        collapsed = ~starved & (z_new < self.collapse_ratio)
        exploded = ~starved & (z_new > self.explode_ratio)
        bad = starved | collapsed | exploded | ~np.isfinite(z_new)
        if bad.any():
            for mask, why in ((starved, f"effective number of points below {self._min_resp:g}"),
                              (collapsed, f"covariance scale fell below {self.collapse_ratio:g}"),
                              (exploded, f"covariance scale exceeded {self.explode_ratio:g}"),
                              (~np.isfinite(z_new), "non-finite scale")):
                if mask.any():
                    self._flag(np.flatnonzero(mask), why, it)
            if self.on_degenerate == "raise":
                raise DegenerateModelError(
                    f"{int(bad.sum())} AGTM component(s) degenerated at iteration {it}: "
                    + "; ".join(sorted({self._reasons[int(k)] for k in np.flatnonzero(bad)})))
            if self.on_degenerate == "stop":
                self._status = "degenerate"
                return False
            self._frozen |= bad
        z_new = np.where(self._frozen, zeta, z_new)
        if self.zeta_bounds is not None:
            z_new = np.clip(z_new, *self.zeta_bounds)
        if not np.all(np.isfinite(W_new)):
            self._flag(range(K), "non-finite weights", it)
            self._status = "degenerate"
            return False
        self.W_, self.means_, self._zeta = W_new, V, z_new
        self._sync()
        return True

    def finish(self, warn: bool = True) -> "AGTM":
        """Close training: publish diagnostics and release the training data."""
        K = self.graph.n_nodes
        self._sync()
        self.loglik_ = np.asarray(self._history)
        self.n_iter_ = max(len(self._history) - 1, 0)
        self.n_samples_ = float(len(self._T) if self._w is None else self._w.sum())
        self.diagnostics_ = {
            "status": self._status, "degenerate": np.asarray(sorted(self._reasons), dtype=int),
            "reasons": self._reasons, "events": self._events, "frozen": self._frozen.copy(),
            "effective_counts": self._G,
        }
        self._T = self._w = self._state = None
        if warn and self._reasons:
            warnings.warn(f"AGTM: {len(self._reasons)} of {K} components flagged as degenerate "
                          f"({'; '.join(sorted(set(self._reasons.values())))}); see diagnostics_",
                          DegenerateModelWarning, stacklevel=3)
        return self

    def fit(self, T, sample_weight=None) -> "AGTM":
        """Fit the model to the noisy points ``T`` (n, d) of this manifold.

        ``sample_weight`` optionally weights the points.
        """
        self.begin(T, sample_weight)
        K = self.graph.n_nodes
        w = self._w
        for it in range(int(self.max_iter) + 1):
            logp, G, S, C = self._estep()
            ll = float(logp.mean() if w is None else (w @ logp) / w.sum())
            hist = self._history
            if not np.isfinite(ll):
                self._flag(range(K), "non-finite log-likelihood", it)
                self._status = "degenerate"
                break
            if hist and ll < hist[-1] - 1e-7 * max(1.0, abs(hist[-1])):
                self._events.append({"iteration": it, "components": [], "reason": "log-likelihood decreased"})
                if self._state is not None:
                    self.W_, self.means_, self._zeta = self._state
                self._status = "loglik_decrease"
                break
            hist.append(ll)
            self._state = (self.W_.copy(), self.means_.copy(), self._zeta.copy())
            if it == self.max_iter:
                break
            if len(hist) > 1 and hist[-1] - hist[-2] < self.tol:
                self._status = "converged"
                break
            if not self._mstep(G, S, C, it):
                break
        return self.finish()

    # ------------------------------------------------------------------ inference
    def _check(self):
        if not hasattr(self, "_Wh"):
            raise RuntimeError("the model is not fitted yet; call fit(T) first")

    def log_density(self, X=None) -> np.ndarray:
        """Log density under the *current* parameters; ``X=None`` uses the training data while
        training is in progress (between ``begin`` and ``finish``)."""
        self._check()
        X = self._T if X is None else as_points(X)
        return _core.agtm_estep(X, self.means_, self._Wh, self._logc, None, False, self._radii)[0]

    def score_samples(self, X) -> np.ndarray:
        """Log density of each row of ``X``."""
        self._check()
        return _core.agtm_estep(as_points(X), self.means_, self._Wh, self._logc, None, False, self._radii)[0]

    def score(self, X) -> float:
        """Mean log-likelihood per point."""
        return float(self.score_samples(X).mean())

    def responsibilities(self, X) -> np.ndarray:
        """Posterior probability of each node for each row of ``X`` -> (n, K)."""
        self._check()
        return _core.agtm_resp(as_points(X), self.means_, self._Wh, self._logc)

    def node_means(self, X, values) -> np.ndarray:
        """Responsibility-weighted mean of ``values`` (n,) or (n, q) at every node (eq. 32)."""
        self._check()
        v = np.asarray(values, float)
        v2 = np.ascontiguousarray(v.reshape(len(v), -1))
        _, G, _, _, Vs = _core.agtm_estep(as_points(X), self.means_, self._Wh, self._logc, v2, True, self._radii)
        out = Vs / np.maximum(G, 1e-300)[:, None]
        return out[:, 0] if v.ndim == 1 else out

    def sample(self, n: int, random_state=None) -> np.ndarray:
        """Draw ``n`` points from the model."""
        self._check()
        rng = rng_from(random_state)
        K, d = self.means_.shape
        k = rng.integers(0, K, n)
        L = np.linalg.cholesky(self.covariances_)
        return self.means_[k] + np.einsum("nab,nb->na", L[k], rng.standard_normal((n, d)))

    def map(self, edge: int, s) -> np.ndarray:
        """Embedding of points at fraction(s) ``s`` along an edge -> (len(s), d)."""
        return (self.W_ @ self.basis(edge, s)).T

    def edge_curvature(self, h: float = 0.01, method: str = "tangent") -> np.ndarray:
        """Directional curvature of the embedding along every edge, evaluated near the first
        node of the edge with five points spaced ``h`` edge lengths apart.

        ``method="tangent"`` uses the change of the unit tangent per unit arc length (eq. 22),
        which is insensitive to the parametrisation speed; ``method="second"`` uses the
        regularised second difference of eq. 23.
        """
        self._check()
        g = self.graph
        s = h * np.arange(1.0, 6.0)
        out = np.empty(g.n_edges)
        if method == "second":
            V = self.means_
            ht = h * np.linalg.norm(V[g.edges[:, 0]] - V[g.edges[:, 1]], axis=1).mean()
        elif method != "tangent":
            raise ValueError("method must be 'tangent' or 'second'")
        for e in range(g.n_edges):
            w = self.map(e, s)
            if method == "second":
                out[e] = np.linalg.norm(w[4] - 2.0 * w[2] + w[0]) / (4.0 * ht * ht)
                continue
            t2, t4 = w[2] - w[0], w[4] - w[2]
            n2, n4 = np.linalg.norm(t2), np.linalg.norm(t4)
            arc = np.linalg.norm(w[2] - w[1]) + np.linalg.norm(w[3] - w[2])
            out[e] = np.linalg.norm(t4 / n4 - t2 / n2) / arc if min(n2, n4, arc) > 0 else 0.0
        return out

    def __repr__(self) -> str:
        fitted = f", iterations={self.n_iter_}, loglik={self.loglik_[-1]:.4g}" if hasattr(self, "zeta_") else ""
        return f"AGTM(dim={self.graph.dim}, nodes={self.graph.n_nodes}, basis={len(self.centres_)}{fitted})"


class MultiManifoldModel:
    """Hierarchical mixture of the AGTMs of all detected manifolds, with an optional uniform
    background component.

    Parameters
    ----------
    models : list of fitted AGTM
    weights : array, optional
        Mixing coefficient of each manifold; default proportional to its number of training points.
    background : float
        Weight of a uniform background over the box ``bounds`` (0 disables it).
    bounds : (lo, hi), optional
        Corners of the background box; required when ``background > 0``.
    """

    def __init__(self, models, weights=None, background: float = 0.0, bounds=None):
        self.models = list(models)
        if not self.models:
            raise ValueError("at least one manifold model is required")
        w = np.asarray([m.n_samples_ for m in self.models] if weights is None else weights, float)
        if len(w) != len(self.models) or (w < 0).any() or w.sum() <= 0:
            raise ValueError("weights must be non-negative, one per model")
        if not 0.0 <= background < 1.0:
            raise ValueError("background must be in [0, 1)")
        self.background = float(background)
        self.weights_ = (1.0 - self.background) * w / w.sum()
        self._log_bg = None
        if self.background > 0:
            if bounds is None:
                raise ValueError("bounds are required for a background component")
            lo, hi = (np.asarray(b, float) for b in bounds)
            self.bounds = (lo, hi)
            self._log_bg = np.log(self.background) - np.log(hi - lo).sum()

    def _parts(self, X) -> np.ndarray:
        X = as_points(X)
        cols = [np.log(w) + m.score_samples(X) if w > 0 else np.full(len(X), -np.inf)
                for m, w in zip(self.models, self.weights_)]
        if self._log_bg is not None:
            lo, hi = self.bounds
            inside = ((X >= lo) & (X <= hi)).all(1)
            cols.append(np.where(inside, self._log_bg, -np.inf))
        return np.column_stack(cols)

    def score_samples(self, X) -> np.ndarray:
        """Log density of each row of ``X`` under the full mixture."""
        return logsumexp(self._parts(X), axis=1)

    def score(self, X) -> float:
        return float(self.score_samples(X).mean())

    def predict_proba(self, X) -> np.ndarray:
        """Posterior membership of each point: one column per manifold (+ background last)."""
        P = self._parts(X)
        return np.exp(P - logsumexp(P, axis=1, keepdims=True))

    def predict(self, X) -> np.ndarray:
        """Most probable manifold of each point (index into ``models``; -1 = background)."""
        lab = np.argmax(self._parts(X), axis=1)
        if self._log_bg is not None:
            lab[lab == len(self.models)] = -1
        return lab

    def sample(self, n: int, random_state=None) -> np.ndarray:
        rng = rng_from(random_state)
        p = np.r_[self.weights_, self.background]
        counts = rng.multinomial(n, p / p.sum())
        parts = [m.sample(c, rng) for m, c in zip(self.models, counts[:-1]) if c]
        if counts[-1]:
            parts.append(rng.uniform(*self.bounds, (counts[-1], len(self.bounds[0]))))
        X = np.vstack(parts)
        return X[rng.permutation(len(X))]

    def __repr__(self) -> str:
        return f"MultiManifoldModel(manifolds={len(self.models)}, background={self.background:g})"
