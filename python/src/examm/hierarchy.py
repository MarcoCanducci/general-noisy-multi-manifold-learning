"""Hierarchical Abstract GTM.

A tree of AGTMs trained top-down with the EM scheme of P. Tino and I. Nabney, "Hierarchical GTM:
constructing localized nonlinear projection manifolds in a principled way", IEEE TPAMI 24(5),
2002, with AGTMs in place of GTMs:

* every model M has a parent-conditional mixing coefficient pi(M | Parent(M));
* the children of a trained parent are fitted by maximising the likelihood restricted to them,
  each point weighted by the (fixed) responsibility of the parent, P(Parent | t);
* siblings compete for points (eq. 39), P(M | t) = P(M | Parent, t) P(Parent | t) (eq. 30), and
  pi(M | Parent) = sum_n P(M | t_n) / sum_n P(Parent | t_n) (eq. 29).

Two additions to the paper, both optional:

* a *residual* leaf per internal node: the parent's own (frozen) model stays as one more sibling
  of its children, so that points of the parent that no child explains remain accounted for;
* a uniform *background* leaf under the root, which absorbs points that no structure explains.

The density of the hierarchy is the mixture of its leaves (eq. 5).
"""
from __future__ import annotations

import warnings

import numpy as np
from scipy.special import logsumexp

from ._util import as_points, rng_from
from .agtm import AGTM, DegenerateModelWarning

__all__ = ["HierarchicalAGTM", "HNode", "build_hierarchy"]


class HNode:
    """A node of the hierarchy: an AGTM (``model``; None for the root) and its place in the tree.

    Attributes
    ----------
    pi : parent-conditional mixing coefficient pi(M | Parent(M))
    pi_residual : coefficient of this node's residual leaf among its children (0 if it has none)
    resp : (n,) responsibility P(M | t_n) of this node for every training point
    history : restricted log-likelihood per EM iteration of the level this node was trained in
    """

    def __init__(self, model=None, parent=None, pi=1.0, label="root"):
        self.model = model
        self.parent = parent
        self.children: list[HNode] = []
        self.pi = float(pi)
        self.pi_residual = 0.0
        self.resp = None
        self.history = None
        self.label = label

    @property
    def level(self) -> int:
        return 1 if self.parent is None else self.parent.level + 1

    @property
    def is_leaf(self) -> bool:
        return not self.children

    @property
    def weight(self) -> float:
        """Unconditional prior: product of the conditional coefficients along the path (eq. 4)."""
        return self.pi if self.parent is None else self.pi * self.parent.weight

    def __repr__(self) -> str:
        m = "virtual root" if self.model is None else f"{self.model.graph.dim}-D, {self.model.graph.n_nodes} nodes"
        return f"HNode({self.label}: {m}, pi={self.pi:.3g}, children={len(self.children)})"


class HierarchicalAGTM:
    """Hierarchy of AGTMs.

    Parameters
    ----------
    background : float
        Initial weight of a uniform background leaf under the root (0 disables it). Its weight
        is re-estimated by EM like any other mixing coefficient.
    residual : bool
        Keep each internal node's own model as a residual leaf among its children.
    residual_init : float
        Initial coefficient of a residual leaf.
    max_iter, tol :
        EM iterations per level and tolerance on the restricted mean log-likelihood.
    eps : float
        Children are trained only on points whose parent responsibility exceeds ``eps``.
    agtm_params : dict
        Keyword arguments for every :class:`AGTM`.
    """

    def __init__(self, *, background: float = 0.01, residual: bool = True, residual_init: float = 0.1,
                 max_iter: int = 30, tol: float = 1e-4, eps: float = 1e-5, agtm_params: dict | None = None):
        if not 0.0 <= background < 1.0:
            raise ValueError("background must be in [0, 1)")
        self.background = background
        self.residual = residual
        self.residual_init = residual_init
        self.max_iter = max_iter
        self.tol = tol
        self.eps = eps
        self.agtm_params = dict(agtm_params or {})

    # ------------------------------------------------------------------ construction
    def set_data(self, X, sample_weight=None, bounds=None) -> "HierarchicalAGTM":
        """Attach the training data. ``bounds`` (lo, hi) is the box of the uniform background;
        default: the bounding box of the data enlarged by 5 %."""
        self.X_ = as_points(X)
        n = len(self.X_)
        self.sample_weight_ = np.ones(n) if sample_weight is None else np.ascontiguousarray(sample_weight, float)
        if self.sample_weight_.shape != (n,) or (self.sample_weight_ < 0).any():
            raise ValueError("sample_weight must hold one non-negative weight per point")
        if bounds is None:
            lo, hi = self.X_.min(0), self.X_.max(0)
            pad = 0.05 * (hi - lo) + 1e-12
            bounds = (lo - pad, hi + pad)
        self.bounds_ = tuple(np.asarray(b, float) for b in bounds)
        self._log_uniform = -float(np.log(self.bounds_[1] - self.bounds_[0]).sum())
        self.root_ = HNode()
        self.root_.resp = np.ones(n)
        self.root_.pi_residual = float(self.background)
        return self

    def _need_data(self):
        if not hasattr(self, "root_"):
            raise RuntimeError("call set_data(X) first")

    def _residual_logp(self, node: HNode, X) -> np.ndarray:
        """Log density of the residual leaf of ``node``: its own model, or the uniform background."""
        if node.model is not None:
            return node.model.score_samples(X)
        lo, hi = self.bounds_
        return np.where(((X >= lo) & (X <= hi)).all(1), self._log_uniform, -np.inf)

    def add_children(self, graphs, parent: HNode | None = None) -> list[HNode]:
        """Fit one AGTM per graph as the children of ``parent`` (default: the root)."""
        self._need_data()
        parent = self.root_ if parent is None else parent
        if parent.children:
            raise ValueError("this node already has children")
        graphs = [g for g in graphs if g.n_nodes >= 2 and g.n_edges >= 1]
        if not graphs:
            return []
        active = np.flatnonzero(parent.resp > self.eps)
        T = np.ascontiguousarray(self.X_[active])
        wpar = parent.resp[active] * self.sample_weight_[active]
        if wpar.sum() <= 0:
            raise ValueError("the parent is responsible for no data")
        has_res = (parent.model is None and self.background > 0) or (parent.model is not None and self.residual)
        pi_res = (parent.pi_residual if parent.model is None else self.residual_init) if has_res else 0.0
        res_logp = self._residual_logp(parent, T) if has_res else None
        models = [AGTM(g, **self.agtm_params).begin(T, sample_weight=wpar) for g in graphs]
        pis = np.full(len(models), (1.0 - pi_res) / len(models))
        history = []
        prev = None
        for it in range(int(self.max_iter) + 1):
            cols = [np.log(max(p, 1e-300)) + m.log_density() for p, m in zip(pis, models)]
            if has_res:
                cols.append(np.log(max(pi_res, 1e-300)) + res_logp)
            A = np.column_stack(cols)
            lse = logsumexp(A, axis=1)
            ll = float(wpar @ np.where(np.isfinite(lse), lse, -1e300) / wpar.sum())
            if history and ll < history[-1] - 1e-7 * max(1.0, abs(history[-1])) and prev is not None:
                for m, st in zip(models, prev[0]):          # roll back to the last good parameters
                    m.W_, m.means_, m._zeta = st
                    m._sync()
                pis, pi_res = prev[1], prev[2]
                break
            history.append(ll)
            if it == self.max_iter or (len(history) > 1 and history[-1] - history[-2] < self.tol):
                break
            prev = ([(m.W_.copy(), m.means_.copy(), m._zeta.copy()) for m in models], pis.copy(), pi_res)
            cond = np.exp(A - lse[:, None])                  # P(M | Parent, t), eq. 39
            cond[~np.isfinite(cond)] = 0.0
            for k, m in enumerate(models):
                w = np.ascontiguousarray(wpar * cond[:, k])  # P(M | t) times the sample weight, eq. 30
                _, G, S, C = m._estep(w)
                m._mstep(G, S, C, it)
                pis[k] = w.sum() / wpar.sum()                # eq. 29
            if has_res:
                pi_res = float((wpar * cond[:, -1]).sum() / wpar.sum())
        # final responsibilities with the final parameters
        cols = [np.log(max(p, 1e-300)) + m.log_density() for p, m in zip(pis, models)]
        if has_res:
            cols.append(np.log(max(pi_res, 1e-300)) + res_logp)
        A = np.column_stack(cols)
        cond = np.exp(A - logsumexp(A, axis=1)[:, None])
        cond[~np.isfinite(cond)] = 0.0
        nodes = []
        for k, m in enumerate(models):
            m._history = list(history)
            m.finish(warn=False)
            node = HNode(m, parent, pis[k], f"{parent.label}.{k}" if parent.model is not None else str(k))
            node.resp = np.zeros(len(self.X_))
            node.resp[active] = parent.resp[active] * cond[:, k]
            node.history = np.asarray(history)
            m.n_samples_ = float((node.resp * self.sample_weight_).sum())
            nodes.append(node)
        parent.children = nodes
        parent.pi_residual = float(pi_res)
        return nodes

    # ------------------------------------------------------------------ structure
    def nodes(self) -> list[HNode]:
        """All nodes with a model, depth first."""
        self._need_data()
        out, stack = [], list(reversed(self.root_.children))
        while stack:
            n = stack.pop()
            out.append(n)
            stack.extend(reversed(n.children))
        return out

    def leaves(self) -> list[tuple]:
        """Leaves of the density: tuples ``(node, kind, weight)`` with kind ``"model"`` (a childless
        node), ``"residual"`` (the frozen model of an internal node) or ``"background"``."""
        self._need_data()
        out = []

        def walk(node, w):
            if node.children:
                for c in node.children:
                    walk(c, w * c.pi)
                if node.pi_residual > 0:
                    out.append((node, "background" if node.model is None else "residual", w * node.pi_residual))
            elif node.model is not None:
                out.append((node, "model", w))

        walk(self.root_, 1.0)
        if not out:
            raise RuntimeError("the hierarchy has no models yet; call add_children first")
        return out

    def leaf_names(self) -> list[str]:
        return ["background" if k == "background" else (n.label + (" (residual)" if k == "residual" else ""))
                for n, k, _ in self.leaves()]

    # ------------------------------------------------------------------ inference
    def _leaf_logp(self, X) -> np.ndarray:
        X = as_points(X)
        cache: dict[int, np.ndarray] = {}
        cols = []
        for node, kind, w in self.leaves():
            if kind == "background":
                lp = self._residual_logp(node, X)
            else:
                if id(node) not in cache:
                    cache[id(node)] = node.model.score_samples(X)
                lp = cache[id(node)]
            cols.append(np.log(max(w, 1e-300)) + lp)
        return np.column_stack(cols)

    def score_samples(self, X) -> np.ndarray:
        """Log density of each row of ``X`` under the hierarchy (mixture of the leaves, eq. 5)."""
        return logsumexp(self._leaf_logp(X), axis=1)

    def score(self, X) -> float:
        return float(self.score_samples(X).mean())

    def leaf_responsibilities(self, X) -> np.ndarray:
        """Posterior of every leaf for every point -> (n, n_leaves); columns as in :meth:`leaves`."""
        A = self._leaf_logp(X)
        R = np.exp(A - logsumexp(A, axis=1, keepdims=True))
        R[~np.isfinite(R)] = 0.0
        return R

    def predict(self, X) -> np.ndarray:
        """Index (into :meth:`leaves`) of the most probable leaf of each point."""
        return np.argmax(self._leaf_logp(X), axis=1)

    def background_probability(self, X) -> np.ndarray:
        """Posterior probability that a point belongs to the uniform background (0 if disabled)."""
        kinds = [k for _, k, _ in self.leaves()]
        if "background" not in kinds:
            return np.zeros(len(as_points(X)))
        return self.leaf_responsibilities(X)[:, kinds.index("background")]

    def residual_probability(self, X) -> np.ndarray:
        """Posterior probability that a point is explained only by the residual of an internal
        node, i.e. it belongs to a structure but to none of its sub-structures."""
        kinds = np.array([k for _, k, _ in self.leaves()])
        return self.leaf_responsibilities(X)[:, kinds == "residual"].sum(1)

    def regional_deviation(self, X, hops: int = 2, min_region: int = 30) -> dict:
        """Regional, across-the-manifold outlier scores from the global model.

        Each point is assigned to its most probable structure among all leaves of the hierarchy
        (the uniform background is left out, so every point gets a structure) and then scored
        with :meth:`AGTM.orthogonal_deviation` on that structure: where it sits on the manifold,
        how far off it lies in the space orthogonal to the local tangent, and how that compares
        with the other points at the same position.

        Returns the arrays of :meth:`AGTM.orthogonal_deviation` plus ``leaf`` (index into
        :meth:`leaves`), ``leaf_name``, ``leaf_probability`` (posterior of that leaf) and
        ``background_probability``.
        """
        X = as_points(X)
        leaves = self.leaves()
        names = np.array(self.leaf_names(), dtype=object)
        R = self.leaf_responsibilities(X)
        structured = np.array([k != "background" for _, k, _ in leaves])
        bgp = R[:, ~structured].sum(1)
        Rs = np.where(structured[None, :], R, -1.0)
        leaf = np.argmax(Rs, axis=1)
        zero = Rs.max(1) <= 0                       # no structure has any posterior: use the densest one
        if zero.any():
            A = self._leaf_logp(X[zero])
            A[:, ~structured] = -np.inf
            leaf[zero] = np.argmax(A, axis=1)
        n = len(X)
        out = {"node": np.full(n, -1, dtype=np.int64), "toward": np.full(n, -1, dtype=np.int64),
               "fraction": np.zeros(n), "distance": np.zeros(n), "tangential": np.zeros(n),
               "d2": np.full(n, np.nan), "dof": np.zeros(n, dtype=np.int64), "p_value": np.full(n, np.nan),
               "regional_rank": np.full(n, np.nan), "region_size": np.zeros(n, dtype=np.int64)}
        for k in np.unique(leaf):
            sel = np.flatnonzero(leaf == k)
            res = leaves[k][0].model.orthogonal_deviation(X[sel], hops=hops, min_region=min_region)
            for key, v in res.items():
                out[key][sel] = v
        out.update(leaf=leaf, leaf_name=names[leaf], leaf_probability=R[np.arange(n), leaf], background_probability=bgp)
        return out

    def regional_outliers(self, X, rank: float = 0.99, p_value: float | None = None, **kw) -> np.ndarray:
        """Mask of points that are among the farthest from their structure at their position:
        ``regional_rank >= rank`` and, if given, a model tail probability below ``p_value``."""
        res = self.regional_deviation(X, **kw)
        flag = res["regional_rank"] >= rank
        if p_value is not None:
            flag &= res["p_value"] < p_value
        return flag

    def sample(self, n: int, random_state=None, include_background: bool = True) -> np.ndarray:
        """Draw ``n`` points from the hierarchy (optionally from its structured leaves only)."""
        rng = rng_from(random_state)
        leaves = [l for l in self.leaves() if include_background or l[1] != "background"]
        p = np.array([w for _, _, w in leaves])
        counts = rng.multinomial(n, p / p.sum())
        parts = []
        for (node, kind, _), c in zip(leaves, counts):
            if c == 0:
                continue
            parts.append(rng.uniform(*self.bounds_, (c, len(self.bounds_[0]))) if kind == "background"
                         else node.model.sample(c, rng))
        X = np.vstack(parts)
        return X[rng.permutation(len(X))]

    def threshold(self, alpha: float = 1e-3, n_samples: int = 100000, random_state=0) -> float:
        """Log-density below which a point is an outlier at level ``alpha``: the ``alpha``-quantile
        of the log density of points generated by the structured leaves themselves (the
        background leaf is left out). The fraction of data flagged is therefore not fixed."""
        if not 0 < alpha < 1:
            raise ValueError("alpha must be in (0, 1)")
        S = self.sample(int(n_samples), random_state, include_background=False)
        return float(np.quantile(self.score_samples(S), alpha))

    def outliers(self, X, alpha: float | None = 1e-3, percentile: float | None = None, **kw) -> np.ndarray:
        """Boolean mask of outliers.

        With ``alpha`` (default) a point is an outlier if its log density is below
        :meth:`threshold`. With ``percentile`` (0-100) the points below that percentile of the
        log density of ``X`` itself are flagged, which always flags that fraction.
        """
        s = self.score_samples(X)
        if percentile is not None:
            return s < np.percentile(s, percentile)
        return s < self.threshold(alpha, **kw)

    def summary(self) -> str:
        self._need_data()
        lines = [f"HierarchicalAGTM: {len(self.X_)} points in {self.X_.shape[1]}-D, {len(self.leaves())} leaves"]

        def walk(node, depth):
            for c in node.children:
                g = c.model.graph
                nd = len(c.model.diagnostics_["degenerate"])
                lines.append("  " * depth + f"- {c.label}: {g.dim}-D, {g.n_nodes} nodes, pi={c.pi:.3f}, "
                             f"weight={c.weight:.3f}, points={c.resp.sum():.0f}"
                             + (f", {nd} degenerate comps" if nd else ""))
                walk(c, depth + 1)
            if node.children and node.pi_residual > 0:
                name = "background (uniform)" if node.model is None else f"residual of {node.label}"
                lines.append("  " * depth + f"- {name}: pi={node.pi_residual:.4f}")

        walk(self.root_, 1)
        return "\n".join(lines)

    def __repr__(self) -> str:
        n = len(self.nodes()) if hasattr(self, "root_") else 0
        return f"HierarchicalAGTM(models={n}, background={self.background:g}, residual={self.residual})"


def build_hierarchy(X, radii, *, preprocess=None, index_preprocess=None, dims=None, min_region: int = 200, min_nodes=None,
                    region_threshold: float = 0.5, sample_weight=None, random_state=None, verbose: bool = False,
                    learner_params: dict | None = None, **hier_params) -> HierarchicalAGTM:
    """Build a hierarchy top-down by scale.

    Level 2 holds the structures detected with ``radii[0]`` on all the data. For each following
    radius, structures are searched with that (finer) radius inside the region of every current
    leaf (the points it is responsible for with probability above ``region_threshold``) and become
    its children.

    Parameters
    ----------
    radii : sequence of float
        Characteristic scale of each level, coarse to fine.
    preprocess : spec or sequence of specs, optional
        Preprocessing for the detection at each level (see :class:`MultiManifoldLearner`); one spec
        is used for all levels. Detection may thin or filter the data, but every model is trained
        on all the points with the responsibilities as weights.
    index_preprocess : spec or sequence of specs, optional
        Separate preprocessing (typically with density thinning) on which the dimensionality index
        is estimated before being carried to all points; see :class:`MultiManifoldLearner`.
    dims : iterable of int or sequence of those, optional
        Intrinsic dimensions searched at each level.
    hier_params :
        Passed to :class:`HierarchicalAGTM` (``background``, ``residual``, ``agtm_params``, ...).

    The learners used for detection are stored in ``hierarchy.detections_``.
    """
    from .pipeline import MultiManifoldLearner

    X = as_points(X)
    radii = [float(r) for r in np.atleast_1d(radii)]
    per_level = lambda v, k: v[k] if isinstance(v, (list, tuple)) and len(v) == len(radii) and not (  # noqa: E731
        len(v) and isinstance(v[0], (int, np.integer))) else v
    H = HierarchicalAGTM(**hier_params).set_data(X, sample_weight)
    H.detections_ = []
    rng = rng_from(random_state)
    frontier = [H.root_]
    for k, r in enumerate(radii):
        nxt = []
        for node in frontier:
            region = np.flatnonzero(node.resp > region_threshold)
            if len(region) < min_region:
                continue
            L = MultiManifoldLearner(r=r, preprocess=per_level(preprocess, k) if preprocess is not None else "saf",
                                     index_preprocess=per_level(index_preprocess, k),
                                     dims=per_level(dims, k), min_nodes=min_nodes, random_state=rng,
                                     **(learner_params or {}))
            try:
                L.set_data(X[region]).run("graphs")
            except ValueError:
                continue
            graphs = [g for gs in L.graphs_.values() for g in gs]
            H.detections_.append({"level": k + 2, "parent": node.label, "r": r, "learner": L, "region": region})
            if not graphs:
                continue
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", DegenerateModelWarning)
                kids = H.add_children(graphs, node)
            nxt.extend(kids)
            if verbose:
                print(f"[examm] level {k + 2}, parent {node.label}: r = {r:g}, {len(region)} points, "
                      f"{len(kids)} structures", flush=True)
        frontier = nxt
        if not frontier:
            break
    return H
