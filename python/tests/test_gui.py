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
    assert a.sv_box.visible and a.w_sv_vals.value and a.w_sv_fvals.value
    # both views are drawn as soon as the filament is selected, each with its own quantities
    assert a.fig_sv_prof.object is not None and a._sv_frames is not None and a.fig_sv_frame.object is not None
    assert set(a._sv_frames["maps"]) - {"density"} == set(a.w_sv_fvals.value)
    before = a.fig_sv_prof.object
    a.w_sv_fvals.value = [a.columns[0]]                    # changing the quantities of one view redraws that view only
    assert set(a._sv_frames["maps"]) - {"density"} == {a.columns[0]} and a.fig_sv_prof.object is before
    a.w_sv_vals.value = [a.columns[1]]
    assert a.fig_sv_prof.object is not before
    a.w_sv_auto.value = False                              # without the box, the buttons decide
    frames = a._sv_frames
    a.w_sv_fvals.value = [a.columns[2]]
    assert a._sv_frames is frames
    a.w_sv_auto.value = True
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


def test_skeleton_and_membership_downloads(app):
    import io
    import zipfile

    z = zipfile.ZipFile(io.BytesIO(app._dl_skeleton().getvalue()))
    assert sorted(z.namelist()) == ["skeleton_edges.csv", "skeleton_nodes.csv"]
    rows = app._dl_membership().getvalue().decode().splitlines()
    assert rows[0].startswith("row,structure,p_on_structure,p_background") and len(rows) == len(app.X) + 1


def test_filter_preview_threshold_radii_and_ball():
    a = App()
    a.w_example.value = "Three filaments (3-D)"
    a._load_example()
    r = a.learner.r
    assert np.isclose(a.w_f_rad.value, 0.5 * r) and np.isclose(a.w_mb_rad.value, 0.7 * r) and np.isclose(a.w_idx_r.value, r)
    a.tabs.active = 2
    # density filter: it works on the diffused points, so the diffusion is the first sub-step
    assert a.pre_tabs._names == ["3a. Diffusion", "3b. Filter"] and a.b_diffuse.visible and not a.b_pre.visible
    a._diffuse_then_filter()
    assert a.pre_tabs.active == 1 and a._scores["points"] is a._diffused[1].diffused
    assert a._preview and a.fig_pre.object is not None and a._ball == (a.w_f_rad.value, "filter radius")
    a.w_f_min.value = 12
    a._on_threshold()
    kept = int(a._threshold()[1].sum())
    assert f"{kept:,}" in a.filter_info.object
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a._run_step_impl(2)
    assert len(a.learner.preprocess_) == kept and not a._preview
    assert [n for n, _ in a.learner.preprocess_.info["steps"]] == ["SAF", "DensityFilter"]      # the paper's order
    a.w_f_mode.value = "core+border"
    border = int(a._threshold()[1].sum())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a._run_step_impl(2)
    assert border > kept and len(a.learner.preprocess_) == border
    a.w_f_mode.value = "core"
    # LAAT: automatic threshold until the scores are computed, then the slider decides
    a.w_filt.value, a.w_diff.value = "laat", "mbms"
    assert a.pre_tabs._names == ["3a. Filter", "3b. Diffusion"] and a.b_pre.visible and not a.b_diffuse.visible
    assert a._scores is None and a._filter_step()[0][1]["threshold"] == "otsu"
    a._filter_scores()
    auto = int(a._threshold()[1].sum())
    v, (thr, keep) = a._scores["values"], a._threshold()
    a.w_la_thr.value = float(np.round(np.log10(np.median(v[keep]) / a._scores["mean"]), 2))   # half of the kept points
    a._on_threshold()
    kept = int(a._threshold()[1].sum())
    assert 0 < kept < auto and isinstance(a._filter_step()[0][1]["threshold"], float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a._run_step_impl(2)
    assert len(a.learner.preprocess_) == kept
    # a radius of its own for the dimensionality index, and the ball of the radius being set
    a.w_idx_r.value = 0.8 * r
    a._on_radius("w_idx_r")
    assert a._ball == (a.w_idx_r.value, "radius of the dimensionality index") and "Scale r" in a.w_layers.value
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a._run_step_impl(3)
    assert np.isclose(a.learner.index_radius, 0.8 * r) and "index_radius=" in a.script()
    a.w_r.value = 1.2 * r                                  # a new scale resets the radii and the scores
    assert np.isclose(a.w_idx_r.value, 1.2 * r) and a._scores is None and a._ball is None


def test_radii_of_graphs_and_density_models_and_ball_placement(app):
    a = app
    r = a.learner.r
    a.w_idx_r.value = 0.9 * r                              # the later steps follow the radius of the index ...
    assert np.isclose(a.w_gr_r.value, 0.9 * r) and np.isclose(a.w_ag_r.value, 0.9 * r)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a._run_step_impl(5)
    L = a.learner
    assert np.isclose(L.radius_of("graphs"), 0.9 * r) and np.isclose(L.radius_of("agtm"), 0.9 * r) and "radii" in L.summary()
    n0 = sum(m.n_points for m in L.manifolds_)
    a.w_ag_r.value = 0.45 * r                              # ... until they are given their own
    a._sync_params()
    assert L.stale_stages == ("agtm",)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a._run_step_impl(5)
    assert sum(m.n_points for m in L.manifolds_) < n0 and "agtm_radius=" in a.script() and "graph_radius=" in a.script()
    a.w_gr_r.value = 1.1 * r
    a._sync_params()
    assert L.stale_stages == ("graphs", "agtm")
    a._on_radius("w_gr_r")
    assert a._ball == (a.w_gr_r.value, "radius of the graphs")
    # the ball is drawn near the boundary of the data, where there are few points
    V = a._view()(a.X)
    c = a._ball_centre(V, a.w_gr_r.value, False)
    lo, hi = V.min(0), V.max(0)
    assert ((c >= lo) & (c <= hi)).all() and (np.minimum(c - lo, hi - c).min() <= 1.3 * a.w_gr_r.value)
    near = (np.linalg.norm(V - c, axis=1) < a.w_gr_r.value).sum()
    assert near <= np.median([(np.linalg.norm(V - p, axis=1) < a.w_gr_r.value).sum() for p in V[::500]])
    a.w_idx_r.value = r
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a._run_step_impl(5)
