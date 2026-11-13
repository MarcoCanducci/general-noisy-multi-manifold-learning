"""Shared helpers."""
from __future__ import annotations

import numpy as np

from . import _core


def as_points(X, name: str = "X") -> np.ndarray:
    """Validate a point cloud and return it as a C-contiguous float64 (N, d) array."""
    A = np.ascontiguousarray(X, dtype=np.float64)
    if A.ndim != 2:
        raise ValueError(f"{name} must be a 2-D array of shape (n_points, n_features), got {A.shape}")
    if A.shape[0] == 0 or A.shape[1] == 0:
        raise ValueError(f"{name} is empty")
    if not np.isfinite(A).all():
        raise ValueError(f"{name} contains NaN or infinite values")
    return A


def check_radius(r, name: str = "r") -> float:
    r = float(r)
    if not np.isfinite(r) or r <= 0:
        raise ValueError(f"{name} must be a positive number, got {r}")
    return r


def build_tree(X: np.ndarray) -> "_core.KDTree":
    return _core.KDTree(X)


def set_num_threads(n: int | None) -> None:
    """Set the number of threads used by the compiled kernels (None, 0 or -1 = all cores)."""
    _core.set_num_threads(0 if n is None or n < 1 else int(n))


def rng_from(random_state) -> np.random.Generator:
    if isinstance(random_state, np.random.Generator):
        return random_state
    return np.random.default_rng(random_state)
