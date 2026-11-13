"""Preprocessing: turn a raw point cloud into a (noisy, diffused) pair of aligned point sets.

The default steps are :class:`SAF` (diffusion) and :class:`DensityFilter`; :class:`LAAT` and
:class:`MBMS`, the filtering and diffusion steps of 1-DREAM, can replace them. Every step
implements the small :class:`Preprocessor` interface, so further denoisers can be plugged in
without touching the rest of the pipeline::

    class MyStep(Preprocessor):
        def apply(self, state, r):
            keep = ...                      # boolean mask over state.noisy
            moved = ...                     # new positions for the kept points
            return state.select(keep).with_diffused(moved)

Steps can be registered by name with :func:`register` or through the ``manydream.preprocessors``
entry-point group of an installed package.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Callable, Sequence

import numpy as np

from ._util import as_points, build_tree, check_radius, rng_from

__all__ = [
    "PreprocessResult", "Preprocessor", "SAF", "DensityFilter", "DensityThin", "LAAT", "MBMS", "Pipeline", "Identity",
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


class DensityThin(Preprocessor):
    """Even out the sampling density by randomly thinning over-dense regions.

    A point with ``n`` neighbours within ``radius`` (itself included) is kept with probability
    ``min(1, target / n)``. Regions with fewer than ``target`` points per ball are left
    untouched; denser regions are reduced to about ``target`` points per ball. Use it as the
    first step when the density of the data spans orders of magnitude, so that one scale ``r``
    can serve the whole data set. The points that are dropped are not lost for good: they can
    be scored or assigned afterwards with the fitted model.

    Parameters
    ----------
    radius : float, optional
        Radius of the ball used to measure density. Default ``radius_factor * r``.
    radius_factor : float
        Used when ``radius`` is None (default 1.0).
    target : float
        Number of points per ball above which a region is thinned.
    random_state : int or Generator
        Seed of the random selection (default 0, so results are reproducible).
    """

    name = "thin"

    def __init__(self, radius: float | None = None, radius_factor: float = 1.0, target: float = 30,
                 random_state=0):
        if target <= 0:
            raise ValueError("target must be positive")
        self.radius = radius
        self.radius_factor = radius_factor
        self.target = target
        self.random_state = random_state

    def apply(self, state, r):
        rad = check_radius(self.radius if self.radius is not None else self.radius_factor * r, "radius")
        counts = build_tree(state.noisy).radius_count(state.noisy, rad)
        prob = np.minimum(1.0, self.target / counts)
        seed = self.random_state
        rng = rng_from(seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed))
        keep = rng.uniform(size=len(counts)) < prob
        out = state.select(keep)
        out.info["thin_counts"] = counts[keep]
        out.info["thin_probability"] = prob[keep]
        out.info["thin_counts_all"] = counts
        return out


def _otsu(v, bins: int = 200) -> float:
    """Value that best splits ``v`` in two groups (largest between-group variance)."""
    h, e = np.histogram(v, bins)
    c = (e[:-1] + e[1:]) / 2
    w, m = np.cumsum(h), np.cumsum(h * c)
    inner = (w > 0) & (w < w[-1])                       # a split must leave points on both sides
    if not inner.any():
        return float(c[0])
    q = np.full(len(c), -np.inf)
    q[inner] = (m[-1] * w[inner] / w[-1] - m[inner]) ** 2 / (w[inner] * (w[-1] - w[inner]))
    k = int(np.argmax(q))
    return float(e[k + 1])                              # upper edge of the last bin of the lower group


class LAAT(Preprocessor):
    """Locally Aligned Ant Technique: keep the points on which random-walking ants leave pheromone.

    After Taghribi et al. (2022), the filtering step of 1-DREAM. Ants jump between neighbouring
    points; a jump from ``i`` to ``j`` has probability proportional to
    ``exp(beta * ((1 - kappa) * F_j + kappa * E_ij))``, where ``F`` is the pheromone on ``j`` and
    ``E_ij`` the alignment of the jump with the principal directions of the neighbourhood of
    ``i`` (weighted by its normalised eigenvalues). Each visit deposits pheromone, which
    evaporates after every epoch. Ants therefore gather on dense, elongated structures, and the
    points with little pheromone are removed. It is an alternative to :class:`DensityFilter`
    that favours aligned structures over merely dense ones; it is most useful for filaments.

    This is a vectorised NumPy implementation of the published rules, not the 1-DREAM C++ code:
    ``F`` and ``E`` are normalised by their largest value in each neighbourhood, and the dynamic
    radius of that code is not included.

    Parameters
    ----------
    radius : float, optional
        Neighbourhood radius of the ants. Default ``radius_factor * r``.
    radius_factor : float
        Used when ``radius`` is None (default 0.7). The radius should be small enough that
        background points have fewer than ``min_neighbours`` neighbours.
    n_ants, n_epochs, n_steps : int
        Number of ants, of epochs, and of jumps of every ant per epoch. Ants are placed at random
        at the start of every epoch.
    kappa : float
        Weight of the alignment against the pheromone, in [0, 1].
    beta : float
        Inverse temperature of the jump probabilities.
    evaporation : float
        Share of the pheromone lost after every epoch.
    min_neighbours : int
        Ants only stand on points with at least this many neighbours.
    threshold : "otsu" or float
        Points with more pheromone are kept. ``"otsu"`` splits the log-pheromone in two groups;
        a number is a multiple of the mean pheromone. The pheromone of every point is returned
        in ``info["laat_pheromone_all"]`` to help choosing it.
    random_state : int or Generator
        Seed of the walks (default 0).
    """

    name = "laat"

    def __init__(self, radius: float | None = None, radius_factor: float = 0.7, n_ants: int = 100,
                 n_epochs: int = 25, n_steps: int = 200, kappa: float = 0.5, beta: float = 5.0,
                 evaporation: float = 0.1, min_neighbours: int = 5, threshold="otsu", random_state=0):
        if not 0 <= kappa <= 1:
            raise ValueError("kappa must lie in [0, 1]")
        if not 0 <= evaporation < 1:
            raise ValueError("evaporation must lie in [0, 1)")
        self.radius = radius
        self.radius_factor = radius_factor
        self.n_ants = n_ants
        self.n_epochs = n_epochs
        self.n_steps = n_steps
        self.kappa = kappa
        self.beta = beta
        self.evaporation = evaporation
        self.min_neighbours = min_neighbours
        self.threshold = threshold
        self.random_state = random_state

    def pheromone(self, X, r: float) -> np.ndarray:
        """Pheromone left on every row of ``X`` (neighbourhood radius ``r``)."""
        X = as_points(X)
        n, d = X.shape
        seed = self.random_state
        rng = rng_from(seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed))
        tree = build_tree(X)
        indptr, idx, dist = tree.radius_neighbors(X, float(r))
        row = np.repeat(np.arange(n), np.diff(indptr))
        other = idx != row
        row, idx, dist = row[other], idx[other], dist[other]
        deg = np.bincount(row, minlength=n)
        indptr = np.r_[0, np.cumsum(deg)]
        ok = deg >= int(self.min_neighbours)
        F = np.full(n, 1e-6)
        if not ok.any() or len(idx) == 0:
            return F
        # alignment of every possible jump with the local principal directions
        _, ev, vec = tree.local_pca(X, float(r), d, 2)
        lam = ev / np.maximum(ev.sum(1, keepdims=True), 1e-300)
        diff = X[idx] - X[row]
        cos = np.abs(np.einsum("ed,ekd->ek", diff, vec[row])) / np.maximum(dist, 1e-300)[:, None]
        E = (cos * lam[row]).sum(1)
        top = np.zeros(n)
        np.maximum.at(top, row, E)
        E /= np.maximum(top[row], 1e-300)
        allowed = ok[row] & ok[idx]                       # ants only stand on well-populated points
        starts = np.flatnonzero(ok)
        deposit = n / float(self.n_ants * self.n_steps)   # about one unit per point and epoch
        first, last = indptr[:-1], indptr[1:] - 1
        for _ in range(int(self.n_epochs)):
            top[:] = 0.0
            np.maximum.at(top, row, F[idx])
            p = np.where(allowed, np.exp(self.beta * ((1 - self.kappa) * F[idx] / top[row] + self.kappa * E - 1.0)), 0.0)
            cum = np.cumsum(p)
            base = np.where(first > 0, cum[np.maximum(first - 1, 0)], 0.0)
            total = np.where(deg > 0, cum[np.maximum(last, 0)] - base, 0.0)
            ants = rng.choice(starts, int(self.n_ants))
            add = np.zeros(n)
            for _ in range(int(self.n_steps)):
                move = total[ants] > 0
                a = ants[move]
                e = np.searchsorted(cum, base[a] + rng.random(len(a)) * total[a], side="right")
                ants[move] = idx[np.clip(e, first[a], last[a])]
                np.add.at(add, ants[move], deposit)
            F = (1 - self.evaporation) * (F + add)
        return F

    def apply(self, state, r):
        rad = check_radius(self.radius if self.radius is not None else self.radius_factor * r, "radius")
        F = self.pheromone(state.noisy, rad)
        if self.threshold == "otsu":
            thr = float(np.exp(_otsu(np.log(F)))) if F.max() > F.min() else float(F.max())
        else:
            thr = float(self.threshold) * float(F.mean())
        keep = F > thr
        out = state.select(keep)
        out.info["laat_pheromone"] = F[keep]
        out.info["laat_pheromone_all"] = F
        out.info["laat_threshold"] = thr
        return out


class MBMS(Preprocessor):
    """Manifold blurring mean shift: move points onto the mean axis (or surface) of a structure.

    After Wang & Carreira-Perpinan (2010), the diffusion step of 1-DREAM. Each iteration moves
    every point towards the Gaussian-weighted mean of its neighbours, but only in the directions
    orthogonal to its local tangent space (the ``dim`` leading principal directions of its
    neighbourhood), so points close in on the structure without sliding along it. It is an
    alternative to :class:`SAF`.

    Few iterations are enough; many iterations slowly pull curved structures towards their
    centre of curvature. This is a vectorised NumPy implementation of the published rule, not
    the 1-DREAM C++ code.

    Parameters
    ----------
    radius : float, optional
        Neighbourhood radius. Default ``radius_factor * r``.
    radius_factor : float
        Used when ``radius`` is None (default 0.7).
    sigma : float, optional
        Scale of the Gaussian weights. Default ``sigma_factor * radius``.
    sigma_factor : float
        Used when ``sigma`` is None (default 0.35).
    dim : int
        Dimension of the tangent space that is left untouched: the intrinsic dimension of the
        structures of interest (1 for filaments).
    n_iter : int
        Number of iterations.
    """

    name = "mbms"

    def __init__(self, radius: float | None = None, radius_factor: float = 0.7, sigma: float | None = None,
                 sigma_factor: float = 0.35, dim: int = 1, n_iter: int = 5):
        if dim < 1:
            raise ValueError("dim must be at least 1")
        self.radius = radius
        self.radius_factor = radius_factor
        self.sigma = sigma
        self.sigma_factor = sigma_factor
        self.dim = dim
        self.n_iter = n_iter

    def apply(self, state, r):
        rad = check_radius(self.radius if self.radius is not None else self.radius_factor * r, "radius")
        sigma = check_radius(self.sigma if self.sigma is not None else self.sigma_factor * rad, "sigma")
        Y = state.diffused
        n, d = Y.shape
        k = int(min(self.dim, d - 1))
        shifts = []
        for _ in range(int(self.n_iter)):
            tree = build_tree(Y)
            indptr, idx, dist = tree.radius_neighbors(Y, rad)          # every point is its own neighbour
            w = np.exp(-dist ** 2 / (2 * sigma ** 2))
            mean = np.add.reduceat(w[:, None] * Y[idx], indptr[:-1], axis=0) / np.add.reduceat(w, indptr[:-1])[:, None]
            shift = mean - Y
            counts, _, vec = tree.local_pca(Y, rad, k, k + 2)
            U = vec[:, :k, :]                                          # (n, k, d) tangent directions
            shift -= np.einsum("nkd,nk->nd", U, np.einsum("nkd,nd->nk", U, shift))
            shift[counts < k + 2] = 0.0
            Y = Y + shift
            shifts.append(float(np.linalg.norm(shift, axis=1).mean()))
        out = state.with_diffused(np.ascontiguousarray(Y))
        out.info["mbms_mean_shift"] = shifts
        out.info["mbms_displacement"] = np.linalg.norm(Y - state.noisy, axis=1)
        return out


class Pipeline(Preprocessor):
    """Apply several preprocessing steps in sequence."""

    name = "pipeline"

    def __init__(self, steps: Sequence[Preprocessor]):
        self.steps = [resolve(s) for s in steps]

    def apply(self, state, r):
        # the rows left after every step are recorded, so that removals can be attributed to steps
        trail = list(state.info.get("steps", []))
        for step in self.steps:
            state = step.apply(state, r)
            if not isinstance(step, Pipeline):
                trail.append((type(step).__name__, state.indices.copy()))
            else:
                trail = list(state.info.get("steps", trail))
            state.info["steps"] = trail
        return state


_REGISTRY: dict[str, Callable[[], Preprocessor]] = {
    "none": Identity,
    "saf": lambda: Pipeline([SAF(), DensityFilter()]),
    "saf-only": SAF,
    "density": DensityFilter,
    "thin": DensityThin,
    "thin+saf": lambda: Pipeline([DensityThin(), SAF(), DensityFilter()]),
    "laat": LAAT,
    "mbms": MBMS,
    "laat+mbms": lambda: Pipeline([LAAT(), MBMS()]),
}


def register(name: str, factory: Callable[[], Preprocessor]) -> None:
    """Register a preprocessing step (or pipeline) under ``name``."""
    _REGISTRY[name] = factory


def _load_entry_points() -> None:
    try:
        from importlib.metadata import entry_points

        for ep in entry_points(group="manydream.preprocessors"):
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
