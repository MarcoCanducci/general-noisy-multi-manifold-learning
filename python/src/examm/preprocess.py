"""Preprocessing: turn a raw point cloud into a (noisy, diffused) pair of aligned point sets.

Every step implements the small :class:`Preprocessor` interface, so alternative denoisers
(e.g. LAAT, MBMS) can be plugged in without touching the rest of the pipeline::

    class MyStep(Preprocessor):
        def apply(self, state, r):
            keep = ...                      # boolean mask over state.noisy
            moved = ...                     # new positions for the kept points
            return state.select(keep).with_diffused(moved)

Steps can be registered by name with :func:`register` or through the ``examm.preprocessors``
entry-point group of an installed package.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Callable, Sequence

import numpy as np

from ._util import as_points, build_tree, check_radius

__all__ = [
    "PreprocessResult", "Preprocessor", "SAF", "DensityFilter", "Pipeline", "Identity",
    "register", "available", "resolve",
]


@dataclass
class PreprocessResult:
    """Row-aligned outcome of preprocessing.

    Attributes
    ----------
    indices : (n,) int array
        Rows of the original input that survived.
    noisy : (n, d) array
        The surviving original points, ``X[indices]``.
    diffused : (n, d) array
        The same points after being moved towards the manifold cores.
    info : dict
        Per-step diagnostics (e.g. per-point scores keyed by step name).
    """

    indices: np.ndarray
    noisy: np.ndarray
    diffused: np.ndarray
    info: dict = field(default_factory=dict)

    @classmethod
    def initial(cls, X) -> "PreprocessResult":
        X = as_points(X)
        return cls(np.arange(X.shape[0]), X, X.copy())

    def select(self, keep) -> "PreprocessResult":
        keep = np.asarray(keep)
        return replace(self, indices=self.indices[keep], noisy=np.ascontiguousarray(self.noisy[keep]),
                       diffused=np.ascontiguousarray(self.diffused[keep]), info=dict(self.info))

    def with_diffused(self, diffused) -> "PreprocessResult":
        diffused = as_points(diffused, "diffused")
        if diffused.shape != self.noisy.shape:
            raise ValueError("diffused points must be row-aligned with the noisy points")
        return replace(self, diffused=diffused, info=dict(self.info))

    def __len__(self) -> int:
        return len(self.indices)


class Preprocessor:
    """Base class of preprocessing steps. Subclasses implement :meth:`apply`."""

    name = "preprocessor"

    def apply(self, state: PreprocessResult, r: float) -> PreprocessResult:  # pragma: no cover
        raise NotImplementedError

    def fit_transform(self, X, r: float) -> PreprocessResult:
        """Run the step on raw data ``X`` with characteristic scale ``r``."""
        return self.apply(PreprocessResult.initial(X), check_radius(r))

    def __repr__(self) -> str:
        args = ", ".join(f"{k}={v!r}" for k, v in vars(self).items() if not k.startswith("_"))
        return f"{type(self).__name__}({args})"


class Identity(Preprocessor):
    """No preprocessing: the data are used both as noisy and as diffused points."""

    name = "none"

    def apply(self, state, r):
        return state


class SAF(Preprocessor):
    """Structure-aware filtering: diffuse points towards high-density regions.

    Each iteration moves every sample to the Gaussian-weighted mean of the original points
    around it, plus a small repulsion between samples that keeps them spread along the
    structure:  ``x <- sum_j w_j q_j / sum_j w_j + mu * repulsion``, ``w = exp(-d^2 / (2 h^2))``.

    Parameters
    ----------
    radius : float, optional
        Kernel scale ``h``. Default ``radius_factor * r``.
    radius_factor : float
        Used when ``radius`` is None (default 0.3). ``h`` should be close to the thickness of
        the noise around the manifolds.
    mu : float
        Repulsion weight in [0, 0.5). Without repulsion surfaces and volumes clump into
        lower-dimensional ridges, which biases the dimensionality index.
    n_iter : int
        Maximum number of iterations.
    tol : float
        Stop early when the mean displacement of an iteration falls below ``tol * h``.
    support : float
        Neighbourhood radius in units of ``h``.
    """

    name = "saf"

    def __init__(self, radius: float | None = None, radius_factor: float = 0.3, mu: float = 0.3,
                 n_iter: int = 10, tol: float = 0.0, support: float = 2.0):
        self.radius = radius
        self.radius_factor = radius_factor
        self.mu = mu
        self.n_iter = n_iter
        self.tol = tol
        self.support = support

    def apply(self, state, r):
        h = check_radius(self.radius if self.radius is not None else self.radius_factor * r, "radius")
        tree = build_tree(state.noisy)
        S = state.diffused
        shifts = []
        for _ in range(int(self.n_iter)):
            new = tree.saf_step(S, h, float(self.mu), float(self.support) * h)
            shifts.append(float(np.linalg.norm(new - S, axis=1).mean()))
            S = new
            if self.tol > 0 and shifts[-1] < self.tol * h:
                break
        out = state.with_diffused(S)
        out.info["saf_mean_shift"] = shifts
        out.info["saf_displacement"] = np.linalg.norm(S - state.noisy, axis=1)
        return out


class DensityFilter(Preprocessor):
    """Remove points lying in sparsely populated regions.

    A point is a *core* point when its ball of radius ``radius`` (default ``radius_factor * r``,
    i.e. ``r / 2``, the SAF kernel scale) holds at least ``min_points`` points both in the noisy
    and in the diffused set. ``min_points`` (the paper's tau) must lie between the typical
    background count and the typical on-manifold count for that radius; the per-point counts
    are returned in ``info`` to help choosing it. With ``mode="core+border"``
    points within ``radius`` of a core point in the diffused set are kept as well; with
    ``mode="core"`` (default) only core points are kept. The rule does not depend on the point order.
    """

    name = "density"

    def __init__(self, radius: float | None = None, min_points: int = 6, mode: str = "core",
                 radius_factor: float = 0.5):
        if mode not in ("core", "core+border"):
            raise ValueError("mode must be 'core' or 'core+border'")
        self.radius = radius
        self.radius_factor = radius_factor
        self.min_points = min_points
        self.mode = mode

    def apply(self, state, r):
        rad = check_radius(self.radius if self.radius is not None else self.radius_factor * r, "radius")
        tn, td = build_tree(state.noisy), build_tree(state.diffused)
        cn = tn.radius_count(state.noisy, rad)
        cd = td.radius_count(state.diffused, rad)
        core = (cn >= self.min_points) & (cd >= self.min_points)
        keep = core
        if self.mode == "core+border":
            keep = td.radius_any(state.diffused, rad, core.astype(np.uint8)).astype(bool)
        out = state.select(keep)
        out.info["density_counts_noisy"] = cn[keep]
        out.info["density_counts_diffused"] = cd[keep]
        out.info["density_core"] = core[keep]
        return out


class Pipeline(Preprocessor):
    """Apply several preprocessing steps in sequence."""

    name = "pipeline"

    def __init__(self, steps: Sequence[Preprocessor]):
        self.steps = [resolve(s) for s in steps]

    def apply(self, state, r):
        for step in self.steps:
            state = step.apply(state, r)
        return state


_REGISTRY: dict[str, Callable[[], Preprocessor]] = {
    "none": Identity,
    "saf": lambda: Pipeline([SAF(), DensityFilter()]),
    "saf-only": SAF,
    "density": DensityFilter,
}


def register(name: str, factory: Callable[[], Preprocessor]) -> None:
    """Register a preprocessing step (or pipeline) under ``name``."""
    _REGISTRY[name] = factory


def _load_entry_points() -> None:
    try:
        from importlib.metadata import entry_points

        for ep in entry_points(group="examm.preprocessors"):
            _REGISTRY.setdefault(ep.name, ep.load())
    except Exception:  # pragma: no cover - plugin discovery must never break the import
        pass


def available() -> list[str]:
    """Names of the registered preprocessing steps."""
    _load_entry_points()
    return sorted(_REGISTRY)


def resolve(spec) -> Preprocessor:
    """Turn ``None``, a name, a list of steps or a step object into a :class:`Preprocessor`."""
    if spec is None:
        return Identity()
    if isinstance(spec, Preprocessor):
        return spec
    if isinstance(spec, str):
        if spec not in _REGISTRY:
            _load_entry_points()
        if spec not in _REGISTRY:
            raise ValueError(f"unknown preprocessing step {spec!r}; available: {available()}")
        return _REGISTRY[spec]()
    if isinstance(spec, (list, tuple)):
        return Pipeline(spec)
    if hasattr(spec, "apply"):
        return spec
    raise TypeError(f"cannot interpret {spec!r} as a preprocessing step")
