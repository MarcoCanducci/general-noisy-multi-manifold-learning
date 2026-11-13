"""One-call estimator chaining preprocessing, dimensionality index, graph extraction and AGTM."""
from __future__ import annotations

import time
import warnings
from dataclasses import dataclass, field

import numpy as np

from ._util import as_points, check_radius, rng_from, set_num_threads
from .agtm import AGTM, DegenerateModelWarning, MultiManifoldModel
from .crawling import ManifoldGraph, multi_manifold_crawl
from .dimension import DimIndexResult, dimensionality_index
from .epsgraph import epsilon_graph
from .preprocess import PreprocessResult, resolve

__all__ = ["MultiManifoldLearner", "Manifold", "STAGES"]

STAGES = ("preprocess", "dimension", "graphs", "agtm")

# which stage each parameter first affects
_PARAM_STAGE = {
    "r": 0, "preprocess": 0,
    "smooth_factor": 1, "dim_min_points": 1, "index_preprocess": 1,
    "dims": 2, "graph": 2, "eta": 2, "beta": 2, "min_points": 2, "min_nodes": 2, "min_residual": 2,
    "angle_threshold": 2, "attach": 2,
    "epsilon": 3, "agtm_params": 3, "background": 3,
    "random_state": 2, "n_jobs": 0, "verbose": 4,
}


@dataclass(eq=False)
class Manifold:
    """One detected manifold.

    Attributes
    ----------
    dim : intrinsic dimension
    graph : ManifoldGraph skeleton (``graph.nodes`` are the crawled positions)
    agtm : fitted AGTM density model (None if fitting was skipped)
    point_indices : rows of the original data attached to this manifold
    """

    dim: int
    graph: ManifoldGraph
    agtm: AGTM | None
    point_indices: np.ndarray
    warnings: list = field(default_factory=list)

    @property
    def n_points(self) -> int:
        return len(self.point_indices)

    def __repr__(self) -> str:
        ll = "" if self.agtm is None else f", loglik={self.agtm.loglik_[-1]:.3g}"
        return (f"Manifold(dim={self.dim}, nodes={self.graph.n_nodes}, edges={self.graph.n_edges}, "
                f"points={self.n_points}{ll})")


class MultiManifoldLearner:
    """Detect and model noisy manifolds of different dimensions in a point cloud.

    Parameters
    ----------
    r : float
        Characteristic scale: radius of the local linearisations used by the dimensionality
        index and by crawling. The single most important parameter.
    preprocess : str, Preprocessor, list or None
        Denoising / diffusion stage. ``"saf"`` (default) runs SAF diffusion followed by the
        density filter; ``None`` uses the data as they are; any object implementing
        :class:`examm.preprocess.Preprocessor` (or a list of them) can be given.
    dims : iterable of int, optional
        Intrinsic dimensions to extract. Default ``1 .. d-1``.
    graph : {"crawl", "epsilon"}
        Manifold crawling (accurate, default) or the faster epsilon-graph alternative.
    eta, beta : float
        Crawling step size and merge tolerance, in units of ``r``.
    epsilon : int
        Spacing in hops of the AGTM basis-function centres.
    smooth_factor : float
        Spatial smoothing radius of the dimensionality index in units of ``r``.
    index_preprocess : preprocessing step, pipeline or name, optional
        If given, the dimensionality index is estimated on the points this preprocessing leaves
        (typically ``Pipeline([DensityThin(...), SAF(), DensityFilter()])``) and carried to every
        modelled point from its nearest neighbour among them. Use it when the density varies so
        much that the smoothed index of the dense core spreads over everything else, while all
        points should still be used to build the graphs. Details are kept in ``index_source_``.
    attach : {"same_dim", "all"}
        Noisy points used to train each AGTM: those within ``r`` of the graph nodes whose
        dimensionality index equals the manifold's dimension, or all of them.
    background : float
        Weight of a uniform background component in the combined model (0 = none).
    agtm_params : dict, optional
        Extra keyword arguments for :class:`examm.AGTM`.
    random_state : int or Generator, optional
    n_jobs : int, optional
        Threads for the compiled kernels (default: all cores).
    verbose : bool

    Attributes (after ``fit``)
    --------------------------
    preprocess_ : PreprocessResult
    dim_result_ : DimIndexResult (for the preprocessed points)
    dim_index_ : (n_samples,) dimensionality index of every input point (0 = removed/undetermined)
    graphs_ : dict dim -> list of ManifoldGraph
    manifolds_ : list of Manifold
    model_ : MultiManifoldModel (None if no manifold was found)
    labels_ : (n_samples,) index of the manifold each input point is attached to (-1 = none)
    timings_ : dict stage -> seconds
    """

    def __init__(self, r: float, *, preprocess="saf", dims=None, graph: str = "crawl", eta: float = 0.75,
                 beta: float = 0.4, epsilon: int = 1, smooth_factor: float = 1.0,
                 dim_min_points: int | None = None, index_preprocess=None, min_points: int = 5, min_nodes: int | None = None,
                 min_residual: int = 20, angle_threshold: float = 2.0 / 3.0, attach: str = "same_dim",
                 background: float = 0.0, agtm_params: dict | None = None, random_state=None,
                 n_jobs: int | None = None, verbose: bool = False):
        self.r = r
        self.preprocess = preprocess
        self.dims = dims
        self.graph = graph
        self.eta = eta
        self.beta = beta
        self.epsilon = epsilon
        self.smooth_factor = smooth_factor
        self.dim_min_points = dim_min_points
        self.index_preprocess = index_preprocess
        self.min_points = min_points
        self.min_nodes = min_nodes
        self.min_residual = min_residual
        self.angle_threshold = angle_threshold
        self.attach = attach
        self.background = background
        self.agtm_params = agtm_params
        self.random_state = random_state
        self.n_jobs = n_jobs
        self.verbose = verbose
        self._done = 0
        self.timings_: dict = {}

    # ------------------------------------------------------------------ parameters
    def get_params(self) -> dict:
        return {k: getattr(self, k) for k in _PARAM_STAGE}

    def set_params(self, **params) -> "MultiManifoldLearner":
        """Change parameters; stages affected by the change are marked stale."""
        for k, v in params.items():
            if k not in _PARAM_STAGE:
                raise ValueError(f"unknown parameter {k!r}; valid: {sorted(_PARAM_STAGE)}")
            setattr(self, k, v)
            self._done = min(self._done, _PARAM_STAGE[k])
        return self

    @property
    def stale_stages(self) -> tuple:
        """Stages that have to be (re)computed."""
        return STAGES[self._done:]

    def _log(self, msg):
        if self.verbose:
            print(f"[examm] {msg}", flush=True)

    def _validate(self):
        check_radius(self.r)
        if self.graph not in ("crawl", "epsilon"):
            raise ValueError("graph must be 'crawl' or 'epsilon'")
        if self.attach not in ("same_dim", "all"):
            raise ValueError("attach must be 'same_dim' or 'all'")

    # ------------------------------------------------------------------ stages
    def run_preprocess(self) -> PreprocessResult:
        t = time.perf_counter()
        self.preprocess_ = resolve(self.preprocess).apply(PreprocessResult.initial(self.X_), float(self.r))
        if len(self.preprocess_) == 0:
            raise ValueError("preprocessing removed every point; decrease min_points or increase r")
        self.timings_["preprocess"] = time.perf_counter() - t
        self._done = 1
        self._log(f"preprocess: kept {len(self.preprocess_)} of {len(self.X_)} points "
                  f"({self.timings_['preprocess']:.2f}s)")
        return self.preprocess_

    def run_dimension(self) -> DimIndexResult:
        if self._done < 1:
            self.run_preprocess()
        t = time.perf_counter()
        if self.index_preprocess is None:
            self.dim_result_ = dimensionality_index(self.preprocess_.diffused, float(self.r),
                                                    smooth_factor=self.smooth_factor,
                                                    min_points=self.dim_min_points)
        else:
            # estimate the index on a separately preprocessed (typically density-thinned) subset,
            # where dense regions cannot outvote the rest, and carry it to every modelled point
            # from its nearest neighbour in that subset (original coordinates)
            from dataclasses import fields
            from scipy.spatial import cKDTree

            sub = resolve(self.index_preprocess).apply(PreprocessResult.initial(self.X_), float(self.r))
            if len(sub) == 0:
                raise ValueError("index_preprocess removed every point")
            res = dimensionality_index(sub.diffused, float(self.r), smooth_factor=self.smooth_factor,
                                       min_points=self.dim_min_points)
            _, nn = cKDTree(sub.noisy).query(self.preprocess_.noisy)
            self.dim_result_ = type(res)(**{f.name: (None if getattr(res, f.name) is None else getattr(res, f.name)[nn])
                                            for f in fields(res)})
            self.index_source_ = {"preprocess": sub, "result": res, "nearest": nn}
        self.dim_index_ = np.zeros(len(self.X_), dtype=int)
        self.dim_index_[self.preprocess_.indices] = self.dim_result_.index
        self.timings_["dimension"] = time.perf_counter() - t
        self._done = 2
        self._log(f"dimension: {self.dim_result_.counts_per_dim()} ({self.timings_['dimension']:.2f}s)")
        return self.dim_result_

    def run_graphs(self) -> dict:
        if self._done < 2:
            self.run_dimension()
        t = time.perf_counter()
        pre, idx = self.preprocess_, self.dim_result_.index
        d = self.X_.shape[1]
        dims = list(range(1, d)) if self.dims is None else [int(j) for j in np.atleast_1d(self.dims)]
        rng = rng_from(self.random_state if not isinstance(self.random_state, np.random.Generator)
                       else self.random_state)
        self.graphs_ = {}
        for j in dims:
            if not 1 <= j <= d:
                raise ValueError(f"dims must lie between 1 and {d}")
            sel = np.flatnonzero(idx == j)
            if len(sel) < max(self.min_residual, 2 * j + 2):
                self.graphs_[j] = []
                continue
            D = np.ascontiguousarray(pre.diffused[sel])
            Q = np.ascontiguousarray(pre.noisy[sel])
            if self.graph == "crawl":
                gs = multi_manifold_crawl(D, self.r, j, eta=self.eta, beta=self.beta, min_points=self.min_points,
                                          min_nodes=self.min_nodes, min_residual=self.min_residual,
                                          random_state=rng)
            else:
                gs = epsilon_graph(D, self.r, j, angle_threshold=self.angle_threshold,
                                   min_points=self.min_points, min_nodes=self.min_nodes, random_state=rng)
            # attach noisy points; indices refer to the preprocessed set
            from ._util import build_tree
            if self.attach == "same_dim":
                tree, back = build_tree(Q), sel
            else:
                tree, back = build_tree(pre.noisy), np.arange(len(pre))
            for g in gs:
                if g.node_index is not None:
                    g.node_index = sel[g.node_index]
                g.points = back[np.flatnonzero(tree.radius_union(g.nodes, float(self.r)))]
            self.graphs_[j] = [g for g in gs if len(g.points) > 0]
        self.timings_["graphs"] = time.perf_counter() - t
        self._done = 3
        self._log("graphs: " + ", ".join(f"{len(v)} of dim {k}" for k, v in self.graphs_.items())
                  + f" ({self.timings_['graphs']:.2f}s)")
        return self.graphs_

    def run_agtm(self) -> list:
        if self._done < 3:
            self.run_graphs()
        t = time.perf_counter()
        pre = self.preprocess_
        params = dict(epsilon=self.epsilon)
        params.update(self.agtm_params or {})
        self.manifolds_ = []
        for j, gs in self.graphs_.items():
            for g in gs:
                model, notes = None, []
                if g.n_edges >= 1 and len(g.points) > g.nodes.shape[1]:
                    with warnings.catch_warnings(record=True) as caught:
                        warnings.simplefilter("always", DegenerateModelWarning)
                        model = AGTM(g, **params).fit(pre.noisy[g.points])
                    notes = [str(w.message) for w in caught]
                self.manifolds_.append(Manifold(j, g, model, pre.indices[g.points], notes))
        self.manifolds_.sort(key=lambda m: (m.dim, -m.n_points))
        fitted = [m for m in self.manifolds_ if m.agtm is not None]
        self.model_ = None
        self.labels_ = np.full(len(self.X_), -1, dtype=int)
        if fitted:
            bounds = (self.X_.min(0), self.X_.max(0)) if self.background > 0 else None
            self.model_ = MultiManifoldModel([m.agtm for m in fitted], background=self.background, bounds=bounds)
            # label attached points by their most probable manifold among those they are attached to
            best = np.full(len(self.X_), -np.inf)
            pos = {id(m): i for i, m in enumerate(fitted)}
            for k, m in enumerate(self.manifolds_):
                if m.agtm is None:
                    continue
                w = self.model_.weights_[pos[id(m)]]
                s = np.log(max(w, 1e-300)) + m.agtm.score_samples(self.X_[m.point_indices])
                upd = s > best[m.point_indices]
                self.labels_[m.point_indices[upd]] = k
                best[m.point_indices[upd]] = s[upd]
        self.timings_["agtm"] = time.perf_counter() - t
        self._done = 4
        self._log(f"agtm: {len(fitted)} models ({self.timings_['agtm']:.2f}s)")
        n_bad = sum(len(m.agtm.diagnostics_["degenerate"]) for m in fitted)
        if n_bad:
            n_all = sum(m.graph.n_nodes for m in fitted)
            warnings.warn(f"{n_bad} of {n_all} mixture components were flagged as degenerate and handled "
                          f"according to on_degenerate; see manifold.agtm.diagnostics_",
                          DegenerateModelWarning, stacklevel=3)
        return self.manifolds_

    # ------------------------------------------------------------------ driver
    def set_data(self, X) -> "MultiManifoldLearner":
        """Attach data without computing anything; use :meth:`run` to compute stage by stage."""
        self.X_ = as_points(X)
        self._done = 0
        self.timings_ = {}
        return self

    def fit(self, X) -> "MultiManifoldLearner":
        """Run the whole pipeline on ``X`` (n_samples, n_features)."""
        self._validate()
        self.set_data(X)
        return self.run(until="agtm")

    def run(self, until: str = "agtm") -> "MultiManifoldLearner":
        """Compute the stale stages up to and including ``until``."""
        if not hasattr(self, "X_"):
            raise RuntimeError("call fit(X) first")
        if until not in STAGES:
            raise ValueError(f"until must be one of {STAGES}")
        self._validate()
        set_num_threads(self.n_jobs)
        steps = (self.run_preprocess, self.run_dimension, self.run_graphs, self.run_agtm)
        for k in range(self._done, STAGES.index(until) + 1):
            steps[k]()
        return self

    def refit(self, from_stage: str | None = None, **params) -> "MultiManifoldLearner":
        """Change parameters and recompute only what depends on them.

        ``from_stage`` forces recomputation from that stage even if no parameter changed.
        """
        self.set_params(**params)
        if from_stage is not None:
            if from_stage not in STAGES:
                raise ValueError(f"from_stage must be one of {STAGES}")
            self._done = min(self._done, STAGES.index(from_stage))
        return self.run()

    # ------------------------------------------------------------------ inference
    def _model(self) -> MultiManifoldModel:
        if getattr(self, "model_", None) is None:
            raise RuntimeError("no fitted manifold model available")
        return self.model_

    def score_samples(self, X) -> np.ndarray:
        """Log density of new points under the combined model."""
        return self._model().score_samples(X)

    def score(self, X) -> float:
        return self._model().score(X)

    def predict_proba(self, X) -> np.ndarray:
        return self._model().predict_proba(X)

    def predict(self, X) -> np.ndarray:
        """Most probable manifold (index into the fitted manifolds; -1 = background)."""
        return self._model().predict(X)

    def summary(self) -> str:
        """Human-readable summary of the fitted pipeline."""
        if not hasattr(self, "X_"):
            return "MultiManifoldLearner (not fitted)"
        n, d = self.X_.shape
        lines = [f"MultiManifoldLearner: {n} points in {d}-D, r = {self.r:g}"]
        if self._done >= 1:
            lines.append(f"  preprocess : kept {len(self.preprocess_)} points ({len(self.preprocess_) / n:.1%})")
        if self._done >= 2:
            c = self.dim_result_.counts_per_dim()
            lines.append("  dimension  : " + ", ".join(f"{v} x {k}-D" if k else f"{v} undetermined"
                                                      for k, v in sorted(c.items())))
        if self._done >= 3:
            lines.append("  graphs     : " + ", ".join(f"{len(v)} of dim {k}" for k, v in self.graphs_.items()))
        if self._done >= 4:
            for k, m in enumerate(self.manifolds_):
                nd = 0 if m.agtm is None else len(m.agtm.diagnostics_["degenerate"])
                flag = f"  [{nd} degenerate component{'s' if nd != 1 else ''}]" if nd else ""
                lines.append(f"  manifold {k:<2d}: {m!r}{flag}")
        if self.stale_stages:
            lines.append(f"  stale      : {', '.join(self.stale_stages)}")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"MultiManifoldLearner(r={self.r!r}, graph={self.graph!r}, preprocess={self.preprocess!r})"
