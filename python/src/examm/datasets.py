"""Synthetic data sets of the paper (Appendix A) and a few analytic test manifolds.

Ground-truth samplers return points on the noiseless manifold. :func:`noisy_sample` thickens
a manifold the way the paper does: a cuboid around the manifold is filled with uniform points
and those within ``radius`` of the manifold are kept.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._util import as_points, build_tree, rng_from

__all__ = [
    "parabolic_arm", "spiral", "s_shape", "hyperbolic_cap", "torus", "mobius", "ball", "torus3",
    "circle", "sphere_bump", "noisy_sample", "background", "d1_ground_truth", "make_d1", "make_three_toroids",
    "make_torus3", "make_filaments", "Dataset",
]

T1 = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]], dtype=float)
T2 = np.array([[0, -1, 0], [0, 0, -1], [1, 0, 0]], dtype=float)


@dataclass
class Dataset:
    """Points with ground truth: ``label`` 0 is background noise, k >= 1 the k-th manifold;
    ``dim`` is the intrinsic dimension of each point's manifold (0 for noise)."""

    X: np.ndarray
    label: np.ndarray
    dim: np.ndarray
    names: list

    def __iter__(self):
        return iter((self.X, self.label))


# ----------------------------------------------------------------- ground truths
def parabolic_arm(n=10000, random_state=None):
    x = rng_from(random_state).uniform(0, 1, n)
    return np.column_stack([x, -2 * x**2 + 2 * x + 1, np.full(n, -0.2)])


def spiral(n=10000, random_state=None):
    s = rng_from(random_state).uniform(0, 1, n)
    theta, rad, z = 2 * np.pi + 8 * np.pi * s, 0.1 + 0.9 * s, 1.0 - 0.9 * s
    return np.column_stack([0.5 * (rad * np.cos(theta) + 1) - 0.1, 0.5 * (rad * np.sin(theta) + 1) + 0.4, z])


def s_shape(n=10000, random_state=None):
    rng = rng_from(random_state)
    th = rng.uniform(-5.0, 3 * np.pi - 5.0, n)
    return np.column_stack([(np.sin(th) + 1) / 5, rng.uniform(1 / 3, 1, n), np.sign(th) * (np.cos(th) + 1) / 3])


def hyperbolic_cap(n=10000, random_state=None):
    rng = rng_from(random_state)
    rad = rng.uniform(0, 0.5, n)
    th = rng.uniform(0, 2 * np.pi, n)
    R = 0.25 - rad**2
    return np.column_stack([R * np.sin(th) + 0.5, R * np.cos(th) - 0.5, 0.4 + 0.4 * rad])


def torus(n=10000, R=0.5, r=0.15, random_state=None):
    rng = rng_from(random_state)
    u, v = rng.uniform(0, 2 * np.pi, n), rng.uniform(0, 2 * np.pi, n)
    return np.column_stack([(R + r * np.cos(v)) * np.cos(u), (R + r * np.cos(v)) * np.sin(u), r * np.sin(v)])


def mobius(n=10000, random_state=None):
    rng = rng_from(random_state)
    u, t = rng.uniform(0, 2 * np.pi, n), rng.uniform(-1, 1, n)
    w = 1 + 0.5 * t * np.cos(u / 2)
    return np.column_stack([w * np.cos(u), w * np.sin(u), 0.5 * t * np.sin(u / 2)])


def ball(n=10000, d=3, radius=1.0, random_state=None):
    rng = rng_from(random_state)
    g = rng.standard_normal((n, d))
    g /= np.linalg.norm(g, axis=1, keepdims=True)
    return radius * g * rng.uniform(0, 1, (n, 1)) ** (1.0 / d)


def torus3(n=20000, random_state=None, grid: int | None = None):
    """3-torus embedded in 4-D (eq. A.6). With ``grid`` a regular grid^3 sampling is returned."""
    if grid:
        a = np.linspace(0, 2 * np.pi, grid, endpoint=False)
        th, ph, ps = (m.ravel() for m in np.meshgrid(a, a, a, indexing="ij"))
    else:
        th, ph, ps = rng_from(random_state).uniform(0, 2 * np.pi, (3, n))
    return np.column_stack([(3 + np.cos(ph)) * np.cos(th), (3 + np.cos(ph)) * np.sin(th),
                            (3 + np.sin(ph)) * np.cos(ps), (3 + np.sin(ph)) * np.sin(ps)])


def circle(n=2000, radius=1.0, d=2, random_state=None):
    th = rng_from(random_state).uniform(0, 2 * np.pi, n)
    X = np.zeros((n, d))
    X[:, 0], X[:, 1] = radius * np.cos(th), radius * np.sin(th)
    return X


def sphere_bump(n=40000, random_state=None):
    """Unit sphere with a Gaussian bump (eq. 24)."""
    rng = rng_from(random_state)
    th, ph = rng.uniform(0, 2 * np.pi, n), rng.uniform(-np.pi / 2, np.pi / 2, n)
    mu, c = np.array([np.pi / 2, -np.pi / 4]), np.array([np.pi / 90, np.pi / 180])
    dens = np.exp(-0.5 * ((th - mu[0]) ** 2 / c[0] + (ph - mu[1]) ** 2 / c[1])) / (2 * np.pi * np.sqrt(c.prod()))
    rad = 1 + 1.5 * dens
    return np.column_stack([rad * np.sin(th) * np.cos(ph), rad * np.sin(th) * np.sin(ph), rad * np.cos(th)])


# ----------------------------------------------------------------- noise models
def noisy_sample(ground_truth, radius=0.04, density=6e4, random_state=None):
    """Uniform points (``density`` per unit volume) lying within ``radius`` of the manifold sample."""
    G = as_points(ground_truth)
    rng = rng_from(random_state)
    lo, hi = G.min(0) - radius, G.max(0) + radius
    U = rng.uniform(lo, hi, (max(int(density * np.prod(hi - lo)), 1), G.shape[1]))
    near = build_tree(G).radius_any(U, radius, np.ones(len(G), np.uint8)).astype(bool)
    return U[near]


def background(X, n, margin=0.2, random_state=None):
    """``n`` uniform points in the bounding box of ``X`` enlarged by ``margin``."""
    X = np.asarray(X)
    return rng_from(random_state).uniform(X.min(0) - margin, X.max(0) + margin, (int(n), X.shape[1]))


def _assemble(parts, dims, names, noise_ratio, margin, rng):
    pts = np.vstack(parts)
    lab = np.concatenate([np.full(len(p), k + 1) for k, p in enumerate(parts)])
    dim = np.concatenate([np.full(len(p), dims[k]) for k, p in enumerate(parts)])
    nb = int(round(noise_ratio * len(pts)))
    if nb > 0:
        B = background(pts, nb, margin, rng)
        pts = np.vstack([pts, B])
        lab = np.concatenate([lab, np.zeros(nb, int)])
        dim = np.concatenate([dim, np.zeros(nb, int)])
    perm = rng.permutation(len(pts))
    return Dataset(np.ascontiguousarray(pts[perm]), lab[perm], dim[perm], ["noise"] + list(names))


def _torus_area_uniform(n, R, r, rng):
    """Area-uniform torus sample by rejection (as in Surf_Toroid.m)."""
    th, ph, w = rng.uniform(0, 2 * np.pi, 2 * n), rng.uniform(0, 2 * np.pi, 2 * n), rng.uniform(0, 1, 2 * n)
    ok = w <= (R + r * np.cos(th)) / (R + r)
    th, ph = th[ok], ph[ok]
    return np.column_stack([(R + r * np.cos(th)) * np.cos(ph), (R + r * np.cos(th)) * np.sin(ph), r * np.sin(th)])


def d1_ground_truth(n=5000, random_state=None) -> list:
    """Noiseless samples of the seven structures of data set D1 (port of MultiManifold_UnDens.m)."""
    rng = rng_from(random_state)
    R1 = T1  # row-vector convention of the MATLAB code: p * R
    x = 0.9 * rng.uniform(size=n)
    parab = np.column_stack([x - 0.1, -2 * x**2 + 2 * x + 1.1, np.full(n, -0.2)])
    t = np.sort(24 * rng.uniform(size=n) + 7)
    rad = np.linspace(0.1, 1.0, n)
    spir = np.column_stack([(rad * np.cos(t) + 1) / 2 - 0.2, (rad * np.sin(t) + 1) / 2, np.linspace(4, 3, n)])
    spir = 0.4 * spir + [0.1, 0.4, -0.9]
    t = 3 * np.pi * rng.uniform(size=n) - 5
    ssh = np.column_stack([(np.sin(t) + 1) / 5, (2 * rng.uniform(size=n) + 1) / 3,
                           (np.sign(t) * (np.cos(t) - 1) + 1) / 3]) @ R1 + [0, 0, 0.5]
    rr = -np.linspace(0, 0.5, n) ** 2 + 0.25
    tb = 2 * np.pi * rng.uniform(size=n)
    cap = 1.2 * np.column_stack([rr * np.sin(tb) + 0.5, rr * np.cos(tb) - 0.5, np.linspace(0.4, 0.6, n)])
    cap += [-0.4, 1.3, 0.3]
    tor = 1.2 * _torus_area_uniform(n, 0.5, 0.15, rng) + [0.2, 0.68, 0.5]
    u, v = rng.uniform(0, 2 * np.pi, n), rng.uniform(-0.5, 0.5, n)
    mob = np.column_stack([(1 + v * np.cos(u / 2)) * np.cos(u), (1 + v * np.cos(u / 2)) * np.sin(u),
                           v * np.sin(u / 2)])
    mob = 0.25 * mob @ R1 @ np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1.0]]) + [0.4, 1.5, -0.2]
    blob = rng.standard_normal((n, 3)) / 20 + [-0.2, -0.5, 0.0]
    return [parab, spir, ssh, cap, tor, mob, blob]


def _from_cloud(truths, dims, names, lo, hi, n_cloud, radius, noise_ratio, rng):
    """Thicken manifolds by keeping the points of one uniform cloud within ``radius`` of them,
    then add background noise as in the MATLAB generator."""
    U = rng.uniform(lo, hi, (int(n_cloud), len(lo)))
    lab = np.zeros(len(U), int)
    for k, G in enumerate(truths):
        near = build_tree(np.ascontiguousarray(G)).radius_any(U, radius, np.ones(len(G), np.uint8)).astype(bool)
        lab[near & (lab == 0)] = k + 1
    pts, lab = U[lab > 0], lab[lab > 0]
    dim = np.asarray(dims)[lab - 1]
    nb = int(round(noise_ratio * len(pts)))
    if nb > 0:
        mx, mn = pts.max(0), pts.min(0)
        mx = mx + 0.2 * (mx - mn)
        mn = mn - 0.2 * (mx - mn)
        B = 1.2 * (mx - mn) * rng.uniform(size=(nb, pts.shape[1])) + mn
        pts, lab, dim = np.vstack([pts, B]), np.r_[lab, np.zeros(nb, int)], np.r_[dim, np.zeros(nb, int)]
    perm = rng.permutation(len(pts))
    return Dataset(np.ascontiguousarray(pts[perm]), lab[perm], dim[perm], ["noise"] + list(names))


def make_d1(scale=1.0, noise_ratio=2.0, radius=0.04, n_truth=5000, random_state=None) -> Dataset:
    """Data set D1 of the paper: two 1-D, four 2-D and one 3-D structures in uniform noise.

    Port of ``DataSet/MultiManifold_UnDens.m``. ``scale`` multiplies the density of the uniform
    cloud the noisy manifolds are cut from (1.0 = the paper's 1e6 points, giving roughly 35 000
    manifold points plus twice as many background points).
    """
    rng = rng_from(random_state)
    names = ["parabolic arm", "spiral", "S-shape", "hyperbolic cap", "torus", "Moebius strip", "blob"]
    return _from_cloud(d1_ground_truth(n_truth, rng), [1, 1, 2, 2, 2, 2, 3], names, [-1.1, -0.9, -1.1],
                       [1.4, 1.8, 1.4], 1e6 * scale, radius, noise_ratio, rng)


def make_three_toroids(scale=1.0, noise_ratio=0.15, radius=0.04, n_truth=5000, random_state=None) -> Dataset:
    """Data set D2 of the paper: three interlocked noisy toroids in uniform noise."""
    rng = rng_from(random_state)
    t = lambda: _torus_area_uniform(n_truth, 0.7, 0.15, rng)  # noqa: E731
    truths = [t() - [0.95, 0, 0], 1.4 * (t() @ T1) + [0, 0, 0.12], t() + [0.95, 0, 0]]
    lo, hi = np.array([-2.0, -1.4, -1.4]), np.array([2.0, 1.4, 1.4])
    return _from_cloud(truths, [2, 2, 2], ["left torus", "middle torus", "right torus"], lo, hi,
                       6e4 * scale * np.prod(hi - lo), radius, noise_ratio, rng)


def make_torus3(n_truth=40000, density=25.0, radius=0.3, noise_ratio=0.0, random_state=None) -> Dataset:
    """Noisy 3-torus embedded in 4-D (paper section 7)."""
    rng = rng_from(random_state)
    pts = noisy_sample(torus3(n_truth, rng), radius, density, rng)
    return _assemble([pts], [3], ["3-torus"], noise_ratio, 0.5, rng)


def make_filaments(density=3e5, noise_ratio=1.0, radius=0.03, random_state=None) -> Dataset:
    """Three 1-D filaments (parabola, spiral, tilted circle) in 3-D uniform noise."""
    rng = rng_from(random_state)
    f1 = noisy_sample(parabolic_arm(8000, rng), radius, density, rng) + [-0.6, 0.0, 0.6]
    f2 = noisy_sample(spiral(8000, rng), radius, density, rng)
    c = circle(8000, 0.45, 3, rng) @ T1.T + [1.0, 0.9, 0.2]
    f3 = noisy_sample(c, radius, density, rng)
    return _assemble([f1, f2, f3], [1, 1, 1], ["parabola", "spiral", "circle"], noise_ratio, 0.15, rng)
