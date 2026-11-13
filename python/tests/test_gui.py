import warnings

import numpy as np
import pytest

pn = pytest.importorskip("panel")
pytest.importorskip("plotly")

from manydream.gui.app import App, suggest_radius  # noqa: E402
from manydream.gui.scene import Scene, render  # noqa: E402


@pytest.fixture(scope="module")
def app():
    pn.extension("plotly")
    a = App()
    a.w_example.value = "Three filaments (3-D)"
    a._load_example()
    return a


def test_walk_through_all_steps(app):
    a = app
    assert a._done() == 2 and a.learner.stale_stages == ("preprocess", "dimension", "graphs", "agtm")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for step in (2, 3, 4, 5):
            a._run_step_impl(step)
            assert a._done() == step + 1 and a.alert.alert_type in ("success", "warning")
    assert len(a.learner.manifolds_) >= 3
    a.w_layers.value = ["Raw data", "Kept points (noisy)", "Diffused points", "Dimensionality index", "Graphs",
                        "Attached points", "AGTM (fitted graph)", "Scale r"]
    a.w_manifold.value = 0
    a.w_quantity.value = "x3"
    a._explore()
    assert a.fig_explore.object is not None
    a.w_pca.value = True
    a.refresh_view()


def test_changing_a_parameter_marks_later_steps_stale(app):
    a = app
    a.w_beta.value = 0.3
    a._sync_params()
    assert a.learner.stale_stages == ("graphs", "agtm") and a._done() == 4
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a._run_step_impl(5)
    assert a.learner.stale_stages == ()


def test_exports_and_generated_script_reproduce_the_session(app, tmp_path):
    a = app
    assert len(a._dl_model().getvalue()) > 1000 and a._dl_summary().getvalue()[:4] == b"\x89PNG"
    rows = a._dl_labels().getvalue().decode().splitlines()
    assert rows[0].startswith("row,kept,dim_index,manifold") and len(rows) == len(a.X) + 1
    ns: dict = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        exec(compile(a.script(), "session.py", "exec"), ns)  # noqa: S102 - our own generated script
    L = ns["learner"]
    assert [(m.dim, m.graph.n_nodes) for m in L.manifolds_] == [(m.dim, m.graph.n_nodes) for m in a.learner.manifolds_]


def test_errors_are_reported_in_the_page(app, tmp_path):
    a = App()
    a._guard(a._run_step_impl, 2)
    assert a.alert.alert_type == "danger" and "load a data set" in a.alert.object
    a.w_path.value = str(tmp_path / "missing.csv")
    a._guard(a._load_path)
    assert a.alert.alert_type == "danger"
    f = tmp_path / "t.csv"
    rng = np.random.default_rng(0)
    T = rng.normal(size=(300, 5))
    T[3, 2] = np.nan
    np.savetxt(f, T, delimiter=",", header="a,b,c,d,e", comments="")
    a.w_path.value = str(f)
    a._guard(a._load_path)
    assert a.alert.alert_type == "success" and a.table.shape == (299, 5) and a.feature_names == ["a", "b", "c"]
    a.w_features.value = ["a"]
    a._guard(a._apply_features)
    assert a.alert.alert_type == "danger"
    a.w_features.value = ["a", "b", "d", "e"]
    a.w_standardise.value = True
    a._guard(a._apply_features)
    assert a.X.shape == (299, 4) and np.allclose(a.X.std(0), 1)


def test_scene_and_radius_helpers(rng):
    X = rng.uniform(size=(5000, 3))
    r = suggest_radius(X)
    import manydream
    assert 20 <= np.median(manydream._core.KDTree(X).radius_count(X[:500], r)) <= 45
    sc = Scene().points("p", X[:100]).points("q", X[:50, :2], np.zeros((50, 3), np.uint8)).lines(
        "l", X[:10], [[0, 1], [1, 2]]).sphere("s", X[0], 0.1)
    assert [l[0] for l in sc.layers] == ["points", "points", "lines", "lines"]
    pane, used = render(sc, "plotly")
    assert used == "plotly" and len(pane.object.data) == 4


def test_laat_mbms_thinned_index_and_viewer_export():
    a = App()
    a.w_example.value = "Three filaments (3-D)"
    a._load_example()
    a.w_filt.value, a.w_diff.value = "laat", "mbms"
    assert a.w_la_ants.visible and a.w_mb_dim.visible and not a.w_saf_h.visible and not a.w_f_min.visible
    a.w_idx_thin.value, a.w_smooth.value = True, 2.0
    assert a.w_idx_rad.visible
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for step in (2, 3, 4, 5):
            a._run_step_impl(step)
    L = a.learner
    assert [n for n, _ in L.preprocess_.info["steps"]] == ["LAAT", "MBMS"] and a.fig_pre.object is not None
    assert hasattr(L, "index_source_") and len(L.manifolds_) >= 1
    script = a.script()
    assert "LAAT(" in script and "MBMS(" in script and "DensityThin(" in script and "index_preprocess=index_preprocess" in script
    ns: dict = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        exec(compile(script, "session.py", "exec"), ns)  # noqa: S102 - our own generated script
    assert [(m.dim, m.graph.n_nodes) for m in ns["learner"].manifolds_] == [(m.dim, m.graph.n_nodes) for m in L.manifolds_]
    html = a._dl_viewer().getvalue().decode()
    assert html.startswith("<!doctype html>") and "Topology (all structures)" in html and "GUI session" in html


def test_filament_views_in_the_explore_step(app):
    a = app
    k = next(i for i, m in enumerate(a.learner.manifolds_) if m.dim == 1 and m.agtm is not None)
    a.w_manifold.value = k
    a._explore()
    assert a.sv_box.visible and a.w_sv_vals.value
    a._sv_profile()
    assert a.fig_sv_prof.object is not None
    a.w_sv_up.value = 1
    a._sv_compute_frames()
    assert a.w_sv_frame.visible and a.w_sv_frame.end == len(a._sv_frames["centers"]) - 1 and a.fig_sv_frame.object is not None
    a.w_sv_frame.value = 0
    a._sv_show_frame()
    a.w_manifold.value = -1
    a._explore()
    assert not a.sv_box.visible and a._sv_frames is None
    a._guard(a._sv_profile)
    assert a.alert.alert_type == "danger"
