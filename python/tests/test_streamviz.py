import numpy as np
import pytest

import manydream
from manydream import streamviz as sv


@pytest.fixture(scope="module")
def stream():
    rng = np.random.default_rng(3)
    n = 6000
    t = rng.uniform(0, 1, n)
    core = np.c_[3 * t, np.sin(2 * t), 0.5 * t ** 2]
    X = core + rng.normal(0, 0.1, (n, 3))
    values = {"radial": np.linalg.norm(X - core, axis=1), "along": t}
    g = manydream.crawl(core, 0.3, 1, seed=int(np.argmin(np.abs(t - 0.5))))
    model = manydream.AGTM(g).fit(X)
    return X, values, model


def test_spine_refines_the_curve_between_its_centres(stream):
    X, values, model = stream
    a, b = sv.spine(model), sv.spine(model, upsample=4)
    K = model.graph.n_nodes
    assert len(a["points"]) == K and len(b["points"]) == K + 4 * (K - 1) and b["is_node"].sum() == K
    assert np.allclose(b["points"][b["is_node"]], a["points"])
    assert np.allclose(np.linalg.norm(b["tangents"], axis=1), 1) and np.all(np.diff(b["arc"]) > 0)


def test_bidimensional_profile_recovers_radial_and_longitudinal_trends(stream):
    X, values, model = stream
    res = sv.bidimensional_profile(model, X, values, r_max=0.3, n_radial=11)
    P = res["profiles"]
    K = model.graph.n_nodes
    assert P["radial"].shape == (K - 1, 10) and set(P) == {"radial", "along", "density"}
    rc = (res["radius"][:-1] + res["radius"][1:]) / 2
    radial = np.nanmean(P["radial"], axis=0)
    assert np.corrcoef(rc, radial)[0, 1] > 0.97                         # grows with the distance from the core
    along = np.nanmean(P["along"], axis=1)
    assert abs(np.corrcoef(np.arange(K - 1), along)[0, 1]) > 0.97       # changes steadily along the curve
    dens = np.nanmean(P["density"], axis=0)
    assert dens[0] > 3 * dens[-1]                                       # the density falls away from the core
    assert sv.plot_bidimensional_profile(res) is not None


def test_comoving_frames_carry_a_consistent_orthonormal_plane(stream, tmp_path):
    X, values, model = stream
    res = sv.comoving_frames(model, X, values, side=0.6, upsample=2, n_pixels=10)
    F = len(res["centers"])
    assert res["pdf"].shape == (F, 11, 11) and res["maps"]["radial"].shape == (F, 11, 11)      # pixels made odd
    u1, u2, t = res["axes"][:, 0], res["axes"][:, 1], res["tangents"]
    for a, b in ((u1, u2), (u1, t), (u2, t)):
        assert np.abs(np.einsum("fd,fd->f", a, b)).max() < 1e-9
    assert np.allclose(np.linalg.norm(u1, axis=1), 1) and np.einsum("fd,fd->f", u1[:-1], u1[1:]).min() > 0.9
    mid = F // 2
    assert np.unravel_index(np.argmax(res["pdf"][mid]), (11, 11)) == (5, 5)                     # the model peaks on the spine
    centre, edge = res["maps"]["radial"][mid, 5, 5], np.nanmean(res["maps"]["radial"][mid, [0, -1], :])
    assert edge > centre
    assert sv.plot_comoving_frame(res, mid, X) is not None
    short = {**res, **{k: res[k][:3] for k in ("centers", "tangents", "axes", "arc", "n_points", "pdf")},
             "maps": {k: v[:3] for k, v in res["maps"].items()}}
    out = sv.animate_comoving_frames(short, tmp_path / "frames.gif", fps=2)
    assert out.exists() and out.stat().st_size > 1000


def test_views_reject_manifolds_that_are_not_curves(rng):
    X = np.c_[rng.uniform(0, 2, (800, 2)), rng.normal(0, 0.01, 800)]
    m = manydream.AGTM(manydream.crawl(X, 0.3, 2, seed=0)).fit(X)
    with pytest.raises(ValueError):
        sv.bidimensional_profile(m, X, r_max=0.2)
