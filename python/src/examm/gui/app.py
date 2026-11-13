"""Step-by-step browser GUI for the multi-manifold pipeline (Panel)."""
from __future__ import annotations

import argparse
import io
import pickle
import tempfile
import traceback
import warnings
from pathlib import Path

import numpy as np

import examm
from examm import datasets, viz
from examm.gui.scene import Scene, pyvista_available, render
from examm.preprocess import LAAT, MBMS, SAF, DensityFilter, DensityThin, Pipeline

STEPS = ["Data", "Features", "Preprocess", "Dimensionality", "Graphs", "AGTM", "Explore & export"]
STAGE_OF_STEP = {2: "preprocess", 3: "dimension", 4: "graphs", 5: "agtm"}
EXAMPLES = {
    "Three interlocked toroids (3-D)": lambda: datasets.make_three_toroids(scale=0.6, random_state=0),
    "Mixed manifolds D1 (3-D)": lambda: datasets.make_d1(random_state=0),
    "Three filaments (3-D)": lambda: datasets.make_filaments(random_state=0),
    "3-torus (4-D)": lambda: datasets.make_torus3(random_state=0),
}
EXAMPLE_R = {"Three interlocked toroids (3-D)": 0.1, "Mixed manifolds D1 (3-D)": 0.1,
             "Three filaments (3-D)": 0.1, "3-torus (4-D)": 1.0}
LAYERS = ["Raw data", "Kept points (noisy)", "Diffused points", "Dimensionality index", "Graphs",
          "Attached points", "AGTM (fitted graph)", "Scale r"]


def suggest_radius(X, target: int = 30, sample: int = 2000, seed: int = 0) -> float:
    """Radius for which a typical point has about ``target`` neighbours."""
    X = np.ascontiguousarray(X, float)
    tree = examm._core.KDTree(X)
    Q = X[np.random.default_rng(seed).choice(len(X), min(sample, len(X)), replace=False)]
    lo, hi = 0.0, float(np.linalg.norm(X.max(0) - X.min(0)))
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if np.median(tree.radius_count(Q, mid)) < target:
            lo = mid
        else:
            hi = mid
    return float(f"{hi:.3g}")


def _colormap(values):
    """Sequential single-hue ramp (light to dark blue) as uint8 RGB."""
    v = np.asarray(values, float)
    lo, hi = np.nanquantile(v, 0.02), np.nanquantile(v, 0.98)
    t = np.clip((v - lo) / (hi - lo + 1e-300), 0, 1)[:, None]
    a, b = np.array([198, 219, 239.0]), np.array([8, 48, 107.0])
    return (a + t * (b - a)).astype(np.uint8)


class App:
    def __init__(self, path: str | None = None, backend: str = "plotly", max_display: int = 40000):
        import panel as pn

        self.pn = pn
        self.backend = "plotly" if backend == "auto" else backend
        self.max_display = max_display
        self.table = None
        self.columns: list[str] = []
        self.learner: examm.MultiManifoldLearner | None = None
        self.X = None
        self.scaling = None
        self.source = ""
        self._build()
        if path:
            self.w_path.value = str(path)
            self._load_path()

    # ------------------------------------------------------------------ layout
    def _build(self):
        pn = self.pn
        W = pn.widgets
        self.status = pn.pane.HTML("", sizing_mode="stretch_width")
        self.alert = pn.pane.Alert("Load a data set to begin.", alert_type="info", sizing_mode="stretch_width")
        self.spinner = pn.indicators.LoadingSpinner(value=False, size=22, color="primary")

        # --- 1 data
        self.w_example = W.Select(name="Built-in example", options=list(EXAMPLES))
        self.b_example = W.Button(name="Load example", button_type="primary")
        self.w_path = W.TextInput(name="File on this machine", placeholder="/path/to/table.csv | .npy | .npz | .h5 | .fits")
        self.b_path = W.Button(name="Load file", button_type="primary")
        self.w_upload = W.FileInput(accept=".csv,.txt,.dat,.tsv,.npy,.npz,.h5,.hdf5,.fits", multiple=False)
        self.data_info = pn.pane.Markdown("")
        self.data_preview = pn.pane.DataFrame(None, max_rows=8, sizing_mode="stretch_width")
        self.b_example.on_click(lambda e: self._guard(self._load_example))
        self.b_path.on_click(lambda e: self._guard(self._load_path))
        self.w_upload.param.watch(lambda e: self._guard(self._load_upload), "value")
        tab_data = pn.Column(
            pn.pane.Markdown("**1. Data** — choose a table of points: rows are points, columns are features."),
            self.w_example, self.b_example, pn.layout.Divider(), self.w_path, self.b_path, pn.layout.Divider(),
            pn.pane.Markdown("Upload a file"), self.w_upload, self.data_info, self.data_preview)

        # --- 2 features
        self.w_features = W.MultiChoice(name="Features to model", options=[], value=[])
        self.w_standardise = W.Checkbox(name="Standardise features (zero mean, unit variance)", value=False)
        self.w_r = W.FloatInput(name="Characteristic scale r", value=0.1, step=0.01, start=1e-12)
        self.b_suggest = W.Button(name="Suggest r", button_type="default")
        self.w_ax = [W.Select(name=f"{a} axis", options=[], width=110) for a in "xyz"]
        self.w_pca = W.Checkbox(name="Show principal components instead", value=False)
        self.w_colour = W.Select(name="Colour raw data by", options=["(none)"], value="(none)")
        self.b_features = W.Button(name="Apply features", button_type="primary")
        self.b_features.on_click(lambda e: self._guard(self._apply_features))
        self.b_suggest.on_click(lambda e: self._guard(self._suggest))
        self.w_r.param.watch(lambda e: self._on_r(), "value")
        for w in self.w_ax + [self.w_pca, self.w_colour]:
            w.param.watch(lambda e: self._guard(self.refresh_view), "value")
        tab_feat = pn.Column(
            pn.pane.Markdown("**2. Features** — pick the columns that span the space in which manifolds are "
                             "searched. All selected features are used in every computation; the axes below "
                             "only control the 3-D view."),
            self.w_features, self.w_standardise, self.b_features, pn.layout.Divider(),
            pn.pane.Markdown("**Scale** — *r* is the radius of the local neighbourhoods used everywhere: it "
                             "should be larger than the noise thickness and smaller than the distance between "
                             "separate structures. Tick *Scale r* under Layers to see it in the view."),
            pn.Row(self.w_r, pn.Column(pn.Spacer(height=18), self.b_suggest)), pn.layout.Divider(),
            pn.pane.Markdown("**View**"), pn.Row(*self.w_ax), self.w_pca, self.w_colour)

        # --- 3 preprocess
        self.w_thin = W.Checkbox(name="Thin over-dense regions first (for data whose density varies a lot)", value=False)
        self.w_thin_rad = W.FloatSlider(name="Thinning radius / r", start=0.2, end=2.0, step=0.05, value=1.0)
        self.w_thin_target = W.IntSlider(name="Thinning: points per ball to keep", start=5, end=200, value=30)
        self.w_filt = W.Select(name="Filter", options={"Density filter": "density", "LAAT ants (1-DREAM)": "laat", "None": "none"})
        self.w_diff = W.Select(name="Diffusion", options={"SAF": "saf", "MBMS (1-DREAM)": "mbms", "None": "none"})
        self.w_saf_h = W.FloatSlider(name="SAF kernel scale h / r", start=0.05, end=1.0, step=0.05, value=0.3)
        self.w_saf_mu = W.FloatSlider(name="SAF repulsion mu", start=0.0, end=0.49, step=0.01, value=0.3)
        self.w_saf_it = W.IntSlider(name="SAF iterations", start=1, end=40, value=10)
        self.w_mb_rad = W.FloatSlider(name="MBMS radius / r", start=0.2, end=1.5, step=0.05, value=0.7)
        self.w_mb_dim = W.IntSlider(name="MBMS: dimension of the structures (1 = filaments)", start=1, end=3, value=1)
        self.w_mb_it = W.IntSlider(name="MBMS iterations", start=1, end=20, value=5)
        self.w_f_rad = W.FloatSlider(name="Filter radius / r", start=0.1, end=2.0, step=0.05, value=0.5)
        self.w_f_min = W.IntSlider(name="Filter: minimum points in the ball", start=1, end=60, value=6)
        self.w_f_mode = W.Select(name="Filter: keep", options={"core points": "core", "core + border points": "core+border"})
        self.w_la_rad = W.FloatSlider(name="LAAT radius of the ants / r", start=0.2, end=1.5, step=0.05, value=0.7)
        self.w_la_ants = W.IntSlider(name="LAAT: number of ants", start=10, end=500, step=10, value=100)
        self.w_la_epochs = W.IntSlider(name="LAAT: epochs", start=1, end=100, value=25)
        self.w_la_steps = W.IntSlider(name="LAAT: steps per epoch", start=10, end=1000, step=10, value=200)
        self.w_la_kappa = W.FloatSlider(name="LAAT: weight of alignment against pheromone", start=0.0, end=1.0, step=0.05, value=0.5)
        self.w_la_thr = W.FloatInput(name="LAAT: pheromone threshold, x mean (0 = automatic)", value=0.0, start=0.0)
        self._pre_groups = {
            "thin": [self.w_thin_rad, self.w_thin_target],
            "saf": [self.w_saf_h, self.w_saf_mu, self.w_saf_it], "mbms": [self.w_mb_rad, self.w_mb_dim, self.w_mb_it],
            "density": [self.w_f_rad, self.w_f_min, self.w_f_mode],
            "laat": [self.w_la_rad, self.w_la_ants, self.w_la_epochs, self.w_la_steps, self.w_la_kappa, self.w_la_thr]}
        for w in (self.w_thin, self.w_filt, self.w_diff):
            w.param.watch(lambda e: self._show_pre_options(), "value")
        self.b_pre = W.Button(name="Run preprocessing", button_type="primary")
        self.b_pre.on_click(lambda e: self._run_step(2))
        self.fig_pre = pn.pane.Matplotlib(None, tight=True, dpi=90, sizing_mode="stretch_width")
        tab_pre = pn.Column(
            pn.pane.Markdown("**3. Preprocess** — a filter removes background points and a diffusion moves the rest "
                             "towards the cores of the structures. The default pair is the density filter with SAF; "
                             "LAAT (ants leaving pheromone on aligned structures) and MBMS (mean shift across the "
                             "local tangent) are the 1-DREAM alternatives, best suited to filaments. LAAT runs before "
                             "the diffusion, the density filter after it. The figure helps choosing the threshold."),
            self.w_thin, self.w_thin_rad, self.w_thin_target, pn.layout.Divider(),
            self.w_filt, self.w_f_rad, self.w_f_min, self.w_f_mode, self.w_la_rad, self.w_la_ants, self.w_la_epochs,
            self.w_la_steps, self.w_la_kappa, self.w_la_thr, pn.layout.Divider(),
            self.w_diff, self.w_saf_h, self.w_saf_mu, self.w_saf_it, self.w_mb_rad, self.w_mb_dim, self.w_mb_it,
            self.b_pre, self.fig_pre)
        self._show_pre_options()

        # --- 4 dimensionality
        self.w_smooth = W.FloatSlider(name="Smoothing radius / r", start=0.0, end=5.0, step=0.1, value=1.0)
        self.w_idx_thin = W.Checkbox(name="Estimate the index on a thinned copy of the data, then carry it to all points",
                                     value=False)
        self.w_idx_rad = W.FloatSlider(name="Thinning radius / r (index only)", start=0.2, end=2.0, step=0.05, value=0.6,
                                       visible=False)
        self.w_idx_target = W.IntSlider(name="Points per ball to keep (index only)", start=5, end=200, value=30, visible=False)
        self.w_idx_weight = W.Checkbox(name="Weight neighbours by the inverse of their density", value=False)

        def _idx_opts(e=None):
            self.w_idx_rad.visible = self.w_idx_target.visible = bool(self.w_idx_thin.value)

        self.w_idx_thin.param.watch(_idx_opts, "value")
        self.b_dim = W.Button(name="Estimate dimensionality", button_type="primary")
        self.b_dim.on_click(lambda e: self._run_step(3))
        self.dim_info = pn.pane.Markdown("")
        self.fig_dim = pn.pane.Matplotlib(None, tight=True, dpi=90, sizing_mode="stretch_width")
        tab_dim = pn.Column(
            pn.pane.Markdown("**4. Dimensionality index** — local PCA in balls of radius *r* gives each diffused "
                             "point the dimension of the structure it lies on. Wider smoothing gives fewer, "
                             "larger structures; narrower smoothing keeps thin features. Where a dense structure "
                             "meets a sparse one, the dense one imposes its dimension on the neighbourhoods around "
                             "the junction: estimating the index on a thinned copy counters that while all points "
                             "are still used to build the graphs. The inverse-density weighting does the same "
                             "without discarding points, but can over-correct when the density varies a lot."),
            self.w_smooth, self.w_idx_thin, self.w_idx_rad, self.w_idx_target, self.w_idx_weight,
            self.b_dim, self.dim_info, self.fig_dim)

        # --- 5 graphs
        self.w_graph = W.Select(name="Method", options={"Manifold crawling": "crawl", "Epsilon-graph (faster)": "epsilon"})
        self.w_dims = W.MultiChoice(name="Dimensions to extract", options=[], value=[])
        self.w_eta = W.FloatSlider(name="Step size eta (x r)", start=0.3, end=1.0, step=0.05, value=0.75)
        self.w_beta = W.FloatSlider(name="Merge tolerance beta (x r)", start=0.1, end=0.7, step=0.05, value=0.4)
        self.w_minnodes = W.IntInput(name="Minimum nodes per graph (0 = automatic)", value=0, start=0)
        self.b_graph = W.Button(name="Extract graphs", button_type="primary")
        self.b_graph.on_click(lambda e: self._run_step(4))
        self.graph_table = pn.pane.DataFrame(None, sizing_mode="stretch_width", index=False)
        tab_graph = pn.Column(
            pn.pane.Markdown("**5. Graphs** — crawling grows a graph skeleton along each structure of the chosen "
                             "dimension; each node has at most 2j edges."),
            self.w_graph, self.w_dims, self.w_eta, self.w_beta, self.w_minnodes, self.b_graph, self.graph_table)

        # --- 6 agtm
        self.w_eps = W.IntSlider(name="Basis spacing epsilon (hops)", start=1, end=6, value=1)
        self.w_iter = W.IntSlider(name="Maximum EM iterations", start=1, end=300, value=100)
        self.w_degen = W.Select(name="If a component degenerates", options={
            "freeze it and continue": "freeze", "stop at the last healthy iteration": "stop", "raise an error": "raise"})
        self.w_zlo = W.FloatInput(name="zeta lower bound (0 = none)", value=0.0, start=0.0, width=180)
        self.w_zhi = W.FloatInput(name="zeta upper bound (0 = none)", value=0.0, start=0.0, width=180)
        self.w_bg = W.FloatSlider(name="Uniform background weight", start=0.0, end=0.9, step=0.05, value=0.0)
        self.b_agtm = W.Button(name="Fit density models", button_type="primary")
        self.b_agtm.on_click(lambda e: self._run_step(5))
        self.agtm_info = pn.pane.Markdown("")
        self.fig_agtm = pn.pane.Matplotlib(None, tight=True, dpi=90, sizing_mode="stretch_width")
        tab_agtm = pn.Column(
            pn.pane.Markdown("**6. Abstract GTM** — a mixture of manifold-aligned Gaussians is fitted on every "
                             "graph, using the noisy points around it."),
            self.w_eps, self.w_iter, self.w_degen, pn.Row(self.w_zlo, self.w_zhi), self.w_bg, self.b_agtm,
            self.agtm_info, self.fig_agtm)

        # --- 7 explore
        self.w_manifold = W.Select(name="Manifold", options={"All": -1}, value=-1)
        self.w_quantity = W.Select(name="Quantity", options=["curvature"], value="curvature")
        self.fig_explore = pn.pane.Matplotlib(None, tight=True, dpi=90, sizing_mode="stretch_width")
        self.dl_model = W.FileDownload(callback=self._dl_model, filename="examm_model.pkl", label="Download fitted model",
                                       button_type="default")
        self.dl_script = W.FileDownload(callback=self._dl_script, filename="examm_session.py",
                                        label="Download Python script", button_type="default")
        self.dl_labels = W.FileDownload(callback=self._dl_labels, filename="examm_labels.csv",
                                        label="Download per-point results (CSV)", button_type="default")
        self.dl_summary = W.FileDownload(callback=self._dl_summary, filename="examm_summary.png",
                                         label="Download summary figure", button_type="default")
        self.dl_viewer = W.FileDownload(callback=self._dl_viewer, filename="examm_viewer.html",
                                        label="Download interactive viewer (HTML)", button_type="primary")
        self.w_manifold.param.watch(lambda e: self._guard(self._explore), "value")
        self.w_quantity.param.watch(lambda e: self._guard(self._explore), "value")
        # the two 1-DREAM views of a filament
        self.w_sv_vals = W.MultiChoice(name="Quantities (columns of the table)", options=[], value=[], max_items=3)
        self.w_sv_rmax = W.FloatInput(name="Profile: largest distance from the core", value=0.1, start=1e-12, width=260)
        self.w_sv_nrad = W.IntSlider(name="Profile: radial bins", start=5, end=40, value=19)
        self.b_sv_prof = W.Button(name="Bi-dimensional profile", button_type="primary")
        self.fig_sv_prof = pn.pane.Matplotlib(None, tight=True, dpi=90, sizing_mode="stretch_width")
        self.w_sv_side = W.FloatInput(name="Frames: side of the plane", value=0.2, start=1e-12, width=260)
        self.w_sv_up = W.IntSlider(name="Frames: points inserted between centres", start=0, end=10, value=5)
        self.w_sv_thick = W.FloatInput(name="Frames: half thickness of the slab (0 = automatic)", value=0.0, start=0.0, width=260)
        self.b_sv_frames = W.Button(name="Co-moving frames", button_type="primary")
        self.w_sv_frame = W.IntSlider(name="Frame", start=0, end=1, value=0, visible=False)
        self.fig_sv_frame = pn.pane.Matplotlib(None, tight=True, dpi=80, sizing_mode="stretch_width")
        self.dl_sv_movie = W.FileDownload(callback=self._dl_frames_movie, filename="comoving_frames.mp4",
                                          label="Download the frames as a movie", button_type="default", visible=False)
        self.sv_info = pn.pane.Markdown("")
        self._sv_frames = None
        self.b_sv_prof.on_click(lambda e: self._guard(self._sv_profile))
        self.b_sv_frames.on_click(lambda e: self._guard(self._sv_compute_frames))
        self.w_sv_frame.param.watch(lambda e: self._guard(self._sv_show_frame), "value_throttled")
        self.sv_box = pn.Column(
            pn.pane.Markdown("**Views of a filament (1-DREAM)** — for a 1-D manifold. The *bi-dimensional profile* "
                             "shows each quantity against the position along the filament and the distance from its "
                             "core. The *co-moving frames* show the quantities, the point density and the fitted "
                             "model on a plane perpendicular to the filament, moved along it."),
            self.w_sv_vals, pn.Row(self.w_sv_rmax, self.w_sv_nrad), self.b_sv_prof, self.fig_sv_prof,
            pn.Row(self.w_sv_side, self.w_sv_up), self.w_sv_thick, self.b_sv_frames, self.sv_info, self.w_sv_frame,
            self.fig_sv_frame, self.dl_sv_movie, visible=False)
        tab_exp = pn.Column(
            pn.pane.Markdown("**7. Explore & export** — select a manifold to isolate it in the view. For 1-D "
                             "manifolds the quantity is shown along the filament, for higher dimensions on the "
                             "unfolded graph."),
            self.w_manifold, self.w_quantity, self.fig_explore, self.sv_box, pn.layout.Divider(),
            pn.pane.Markdown("The **interactive viewer** is a single HTML page that opens in any browser, without "
                             "this app: the data and the structures in 3-D, each structure's unfolded graph, all "
                             "structures together, outliers by percentile or threshold for several criteria, and "
                             "the columns you did not model as colours. Building it takes up to a few minutes."),
            self.dl_viewer, self.dl_script, self.dl_model, self.dl_labels, self.dl_summary)

        self.tabs = pn.Tabs(*[(f"{k + 1}. {n}", t) for k, (n, t) in enumerate(zip(
            STEPS, [tab_data, tab_feat, tab_pre, tab_dim, tab_graph, tab_agtm, tab_exp]))],
            tabs_location="left", dynamic=False, width=560)
        self.tabs.param.watch(lambda e: self._guard(self._on_tab), "active")

        self.w_layers = W.CheckBoxGroup(name="Layers", options=LAYERS, value=["Raw data"], inline=True)
        self.w_layers.param.watch(lambda e: self._guard(self.refresh_view), "value")
        self.view_holder = pn.Column(pn.pane.Markdown("*The 3-D view appears here once data are loaded.*"),
                                     sizing_mode="stretch_width", min_height=640)
        self.summary = pn.pane.Markdown("", sizing_mode="stretch_width")
        self.backend_note = pn.pane.Markdown("", styles={"color": "#52514e", "font-size": "11px"})
        self.w_az = W.IntSlider(name="Azimuth", start=-180, end=180, step=5, value=35, width=220)
        self.w_el = W.IntSlider(name="Elevation", start=-85, end=85, step=5, value=25, width=220)
        self.w_zoom = W.FloatSlider(name="Zoom", start=0.5, end=6.0, step=0.1, value=1.0, width=220)
        for w in (self.w_az, self.w_el, self.w_zoom):
            w.param.watch(lambda e: self._guard(self.refresh_view), "value_throttled")
        self.cam_row = pn.Row(self.w_az, self.w_el, self.w_zoom, visible=self.backend == "pyvista")
        right = pn.Column(pn.Row(pn.pane.Markdown("**Layers**"), self.w_layers), self.cam_row, self.view_holder,
                          self.backend_note, self.summary, sizing_mode="stretch_width")
        self.root = pn.Column(
            pn.Row(pn.pane.Markdown("## examm — noisy multi-manifold learning"), pn.Spacer(), self.spinner),
            self.status, self.alert, pn.Row(self.tabs, right, sizing_mode="stretch_width"),
            sizing_mode="stretch_width")
        self._update_status()

    def layout(self):
        return self.root

    # ------------------------------------------------------------------ helpers
    def _say(self, text, kind="info"):
        self.alert.object = text
        self.alert.alert_type = kind

    def _guard(self, fn, *a):
        """Run a callback, reporting errors in the page instead of the server log only."""
        self.spinner.value = True
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                out = fn(*a)
            notes = [str(w.message) for w in caught if issubclass(w.category, examm.DegenerateModelWarning)]
            if notes and self.alert.alert_type != "danger":
                self._say(self.alert.object + "<br>Note: " + " | ".join(dict.fromkeys(notes)), "warning")
            return out
        except Exception as exc:  # noqa: BLE001 - shown to the user
            traceback.print_exc()
            self._say(f"**{type(exc).__name__}**: {exc}", "danger")
        finally:
            self.spinner.value = False

    def _done(self) -> int:
        """Number of completed steps (0..7)."""
        if self.table is None:
            return 0
        if self.learner is None:
            return 1
        return 2 + self.learner._done

    def _update_status(self):
        done = self._done()
        chips = []
        for k, name in enumerate(STEPS):
            if k < done or (k == 6 and done >= 6):
                bg, fg, mark = "#e6f4ea", "#0b5a0b", "✓ "
            elif k == done:
                bg, fg, mark = "#e8f0fc", "#17457f", "▶ "
            else:
                bg, fg, mark = "#f1f0ec", "#52514e", ""
            chips.append(f'<span style="background:{bg};color:{fg};padding:3px 10px;border-radius:12px;'
                         f'font-size:12px;margin-right:6px;white-space:nowrap">{mark}{k + 1}. {name}</span>')
        self.status.object = '<div style="line-height:2">' + "".join(chips) + "</div>"
        self.summary.object = ("```\n" + self.learner.summary() + "\n```") if self.learner is not None else ""

    # ------------------------------------------------------------------ data
    def _set_table(self, A, names, source, r=None):
        A = np.asarray(A, float)
        ok = np.isfinite(A).all(1)
        dropped = int((~ok).sum())
        A = A[ok]
        if len(A) < 10:
            raise ValueError("the table has fewer than 10 complete rows")
        self.table, self.columns, self.source = A, list(names), source
        self.learner = None
        import pandas as pd

        self.data_preview.object = pd.DataFrame(A[:8], columns=self.columns).round(4)
        self.data_info.object = (f"**{source}**: {len(A):,} rows, {A.shape[1]} columns"
                                 + (f" ({dropped:,} rows with missing values dropped)" if dropped else ""))
        self.w_features.options = self.columns
        self.w_features.value = self.columns[:min(3, len(self.columns))] if r is None else self.columns
        self.w_colour.options = ["(none)"] + self.columns
        self.w_colour.value = "(none)"
        self.w_quantity.options = ["curvature"] + self.columns
        self._apply_features(r=r)
        self._say(f"Loaded {len(A):,} points. Check the features and the scale *r* in step 2, then run the "
                  "steps in order.", "success")

    def _load_example(self):
        name = self.w_example.value
        ds = EXAMPLES[name]()
        self._set_table(ds.X, [f"x{k + 1}" for k in range(ds.X.shape[1])], name, r=EXAMPLE_R[name])

    def _load_path(self):
        path = Path(self.w_path.value.strip()).expanduser()
        if not path.is_file():
            raise FileNotFoundError(f"no such file: {path}")
        A, names = examm.io.load_table(path)
        self._set_table(A, names, path.name)

    def _load_upload(self):
        if not self.w_upload.value:
            return
        suffix = Path(self.w_upload.filename or "upload.csv").suffix
        with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as fh:
            fh.write(self.w_upload.value)
        try:
            A, names = examm.io.load_table(fh.name)
        finally:
            Path(fh.name).unlink(missing_ok=True)
        self._set_table(A, names, self.w_upload.filename or "upload")

    # ------------------------------------------------------------------ features
    def _apply_features(self, r=None):
        if self.table is None:
            raise ValueError("load a data set first")
        feats = list(self.w_features.value)
        if len(feats) < 2:
            raise ValueError("select at least two features")
        cols = [self.columns.index(f) for f in feats]
        X = self.table[:, cols].copy()
        self.scaling = None
        if self.w_standardise.value:
            mu, sd = X.mean(0), X.std(0)
            sd[sd == 0] = 1.0
            X = (X - mu) / sd
            self.scaling = (mu, sd)
        self.X = np.ascontiguousarray(X)
        self.feature_names = feats
        if r is None:
            r = suggest_radius(self.X)
        self.learner = examm.MultiManifoldLearner(r=r, random_state=0).set_data(self.X)
        self._r_silent = True
        self.w_r.value = r
        self.w_r.step = float(f"{r / 10:.1g}")
        self._r_silent = False
        for k, w in enumerate(self.w_ax):
            w.options = feats + (["(none)"] if k == 2 else [])
            w.value = feats[k] if k < len(feats) else "(none)"
        d = len(feats)
        self.w_dims.options = list(range(1, d + 1))
        self.w_dims.value = list(range(1, d))
        self.w_layers.value = ["Raw data"]
        self.w_manifold.options = {"All": -1}
        self.refresh_view()
        self._update_status()

    def _suggest(self):
        if self.X is None:
            raise ValueError("load a data set first")
        self.w_r.value = suggest_radius(self.X)
        self._say(f"Suggested r = {self.w_r.value:g}: a typical point has about 30 neighbours within it. "
                  "Treat it as a starting point.", "info")

    def _on_r(self):
        if self.learner is None or getattr(self, "_r_silent", False):
            return
        if self.w_r.value and self.w_r.value > 0 and self.w_r.value != self.learner.r:
            self.learner.set_params(r=float(self.w_r.value))
            self._update_status()
            self._guard(self.refresh_view)

    # ------------------------------------------------------------------ pipeline
    def _show_pre_options(self):
        on = {"thin": bool(self.w_thin.value), "saf": self.w_diff.value == "saf", "mbms": self.w_diff.value == "mbms",
              "density": self.w_filt.value == "density", "laat": self.w_filt.value == "laat"}
        for key, widgets in self._pre_groups.items():
            for w in widgets:
                w.visible = on[key]

    def _steps(self, thin):
        """Preprocessing steps as (object, code) pairs; ``thin`` is None or (radius factor, target)."""
        code = lambda cls, kw: cls.__name__ + "(" + ", ".join(f"{k}={v!r}" for k, v in kw.items()) + ")"  # noqa: E731
        out = []

        def add(cls, **kw):
            out.append((cls(**kw), code(cls, kw)))

        if thin is not None:
            add(DensityThin, radius_factor=thin[0], target=thin[1], random_state=0)
        if self.w_filt.value == "laat":
            add(LAAT, radius_factor=self.w_la_rad.value, n_ants=self.w_la_ants.value, n_epochs=self.w_la_epochs.value,
                n_steps=self.w_la_steps.value, kappa=self.w_la_kappa.value,
                threshold=self.w_la_thr.value if self.w_la_thr.value > 0 else "otsu")
        if self.w_diff.value == "saf":
            add(SAF, radius_factor=self.w_saf_h.value, mu=self.w_saf_mu.value, n_iter=self.w_saf_it.value)
        elif self.w_diff.value == "mbms":
            add(MBMS, radius_factor=self.w_mb_rad.value, dim=self.w_mb_dim.value, n_iter=self.w_mb_it.value)
        if self.w_filt.value == "density":
            add(DensityFilter, radius_factor=self.w_f_rad.value, min_points=self.w_f_min.value, mode=self.w_f_mode.value)
        return out

    @staticmethod
    def _as_spec(steps):
        return None if not steps else (steps[0][0] if len(steps) == 1 else Pipeline([s for s, _ in steps]))

    @staticmethod
    def _as_code(steps):
        return "None" if not steps else (steps[0][1] if len(steps) == 1 else "Pipeline([" + ", ".join(c for _, c in steps) + "])")

    def _main_steps(self):
        return self._steps((self.w_thin_rad.value, self.w_thin_target.value) if self.w_thin.value else None)

    def _index_steps(self):
        return self._steps((self.w_idx_rad.value, self.w_idx_target.value)) if self.w_idx_thin.value else []

    def _preprocess_spec(self):
        return self._as_spec(self._main_steps())

    def _params(self) -> dict:
        zb = None
        if self.w_zlo.value > 0 or self.w_zhi.value > 0:
            zb = (self.w_zlo.value if self.w_zlo.value > 0 else 0.0, self.w_zhi.value if self.w_zhi.value > 0 else np.inf)
        return dict(
            r=float(self.w_r.value), preprocess=self._preprocess_spec(), smooth_factor=self.w_smooth.value,
            index_preprocess=self._as_spec(self._index_steps()), index_density_weight=bool(self.w_idx_weight.value),
            graph=self.w_graph.value, dims=sorted(int(j) for j in self.w_dims.value), eta=self.w_eta.value,
            beta=self.w_beta.value, min_nodes=self.w_minnodes.value or None, epsilon=self.w_eps.value,
            background=self.w_bg.value,
            agtm_params=dict(max_iter=self.w_iter.value, on_degenerate=self.w_degen.value, zeta_bounds=zb))

    def _sync_params(self):
        """Push widget values to the learner; only real changes mark stages stale."""
        L = self.learner
        for k, v in self._params().items():
            if repr(getattr(L, k)) != repr(v):
                L.set_params(**{k: v})

    def _run_step(self, step: int):
        self._guard(self._run_step_impl, step)

    def _run_step_impl(self, step: int):
        if self.learner is None:
            raise ValueError("load a data set and apply the features first")
        if step == 4 and not self.w_dims.value:
            raise ValueError("select at least one dimension to extract")
        if self.w_beta.value >= self.w_eta.value:
            raise ValueError("the merge tolerance beta must be smaller than the step size eta")
        self._sync_params()
        L = self.learner
        self._say(f"Running {STEPS[step].lower()}…", "info")
        L.run(until=STAGE_OF_STEP[step])
        t = L.timings_.get(STAGE_OF_STEP[step], 0.0)
        if step == 2:
            n = len(L.preprocess_)
            self._say(f"Preprocessing kept **{n:,}** of {len(self.X):,} points ({n / len(self.X):.1%}) in {t:.1f} s.",
                      "success")
            self.w_layers.value = ["Raw data", "Diffused points"]
        elif step == 3:
            c = L.dim_result_.counts_per_dim()
            self._say("Dimensionality index: " + ", ".join(
                f"**{v:,}** points {k}-D" if k else f"{v:,} undetermined" for k, v in sorted(c.items()))
                + f" ({t:.1f} s).", "success")
            self.w_layers.value = ["Dimensionality index"]
        elif step == 4:
            n = sum(len(v) for v in L.graphs_.values())
            self._say(f"Found **{n}** manifold graph{'s' if n != 1 else ''}: " + ", ".join(
                f"{len(v)} of dimension {k}" for k, v in L.graphs_.items()) + f" ({t:.1f} s).",
                "success" if n else "warning")
            self.w_layers.value = ["Kept points (noisy)", "Graphs"]
        elif step == 5:
            n = len([m for m in L.manifolds_ if m.agtm is not None])
            self._say(f"Fitted **{n}** density model{'s' if n != 1 else ''} in {t:.1f} s.", "success" if n else "warning")
            self.w_layers.value = ["Attached points", "AGTM (fitted graph)"]
            self.w_manifold.options = {"All": -1, **{f"{k}: {m.dim}-D, {m.graph.n_nodes} nodes, {m.n_points:,} points": k
                                                     for k, m in enumerate(L.manifolds_)}}
            self.w_manifold.value = -1
        self._diagnostics(step)
        self.refresh_view()
        self._update_status()

    def _fig(self, pane, fig):
        import matplotlib.pyplot as plt

        pane.object = fig
        plt.close(fig)

    def _diagnostics(self, step: int):
        import pandas as pd

        L = self.learner
        if step == 2:
            if self.w_filt.value == "density":
                self._fig(self.fig_pre, viz.plot_density_counts(self.X, self.w_f_rad.value * L.r, self.w_f_min.value))
            elif self.w_filt.value == "laat" and "laat_pheromone_all" in L.preprocess_.info:
                import matplotlib.pyplot as plt

                F, thr = L.preprocess_.info["laat_pheromone_all"], L.preprocess_.info["laat_threshold"]
                fig, ax = plt.subplots(figsize=(5.5, 3.6))
                ax.hist(np.log10(F), bins=80, color=viz.CATEGORICAL[0], edgecolor="white", linewidth=0.2)
                ax.axvline(np.log10(thr), color="#0b0b0b", linewidth=1, linestyle="--")
                ax.set(xlabel="log10 pheromone on a point", ylabel="points", title=f"kept: pheromone above {thr:.3g}")
                ax.spines[["top", "right"]].set_visible(False)
                self._fig(self.fig_pre, fig)
            else:
                self.fig_pre.object = None
        elif step == 3:
            c = L.dim_result_.counts_per_dim()
            self.dim_info.object = " | ".join(f"**{k}-D**: {v:,}" if k else f"undetermined: {v:,}" for k, v in sorted(c.items()))
            self._fig(self.fig_dim, viz.plot_simplex(L.dim_result_, figsize=(7, 3.4)))
        elif step == 4:
            rows = [dict(dim=j, nodes=g.n_nodes, edges=g.n_edges, points=len(g.points), length=round(g.total_length(), 3))
                    for j, gs in L.graphs_.items() for g in gs]
            self.graph_table.object = pd.DataFrame(rows) if rows else pd.DataFrame({"result": ["no graphs found"]})
        elif step == 5:
            fitted = [m for m in L.manifolds_ if m.agtm is not None]
            if fitted:
                self._fig(self.fig_agtm, viz.plot_loglik(L, figsize=(6, 3.4)))
            lines = []
            for k, m in enumerate(L.manifolds_):
                if m.agtm is None:
                    continue
                dg = m.agtm.diagnostics_
                lines.append(f"- manifold {k}: {m.agtm.n_iter_} iterations, {dg['status']}, "
                             f"{len(dg['degenerate'])} of {m.graph.n_nodes} components degenerate")
            self.agtm_info.object = "\n".join(lines)

    def _on_tab(self):
        if self.tabs.active == 6:
            self._explore()

    def _explore(self):
        L = self.learner
        if L is None or L._done < 4 or not L.manifolds_:
            self.fig_explore.object = None
            return
        k = self.w_manifold.value
        self.refresh_view()
        is_curve = k is not None and k >= 0 and L.manifolds_[k].dim == 1 and L.manifolds_[k].agtm is not None
        if self.sv_box.visible != is_curve or getattr(self, "_sv_for", None) != k:
            self.sv_box.visible, self._sv_for, self._sv_frames = is_curve, k, None
            self.fig_sv_prof.object = self.fig_sv_frame.object = None
            self.w_sv_frame.visible = self.dl_sv_movie.visible = False
            self.sv_info.object = ""
            if is_curve:
                passive = [c for c in self.columns if c not in self.feature_names]
                self.w_sv_vals.options = list(self.columns)
                self.w_sv_vals.value = (passive or list(self.feature_names))[:2]
                self.w_sv_rmax.value, self.w_sv_side.value = float(L.r), 2.0 * float(L.r)
        if k is None or k < 0:
            self._fig(self.fig_explore, viz.plot_loglik(L, figsize=(6, 3.4)))
            return
        m = L.manifolds_[k]
        if m.agtm is None:
            self.fig_explore.object = None
            return
        q = self.w_quantity.value
        vals = None if q == "curvature" else self.table[:, self.columns.index(q)]
        if m.dim == 1:
            fig = viz.plot_profile(m, self.X if vals is not None else None, vals, label=q, figsize=(6, 3.2))
        elif vals is None:
            fig = viz.plot_unfolded(m, m.agtm.edge_curvature(), label="curvature", figsize=(5.5, 4.6))
        else:
            fig = viz.plot_unfolded(m, node_values=m.agtm.node_means(self.X, vals), label=q, figsize=(5.5, 4.6))
        self._fig(self.fig_explore, fig)

    # ------------------------------------------------------------------ views of a filament
    def _sv_target(self):
        L = self._need_fit()
        k = self.w_manifold.value
        if k is None or k < 0 or L.manifolds_[k].dim != 1 or L.manifolds_[k].agtm is None:
            raise ValueError("select a fitted 1-D manifold first")
        return L.manifolds_[k], {c: self.table[:, self.columns.index(c)] for c in self.w_sv_vals.value}

    def _sv_profile(self):
        from examm import streamviz as sv

        m, vals = self._sv_target()
        res = sv.bidimensional_profile(m, self.X, vals, r_max=float(self.w_sv_rmax.value), n_radial=int(self.w_sv_nrad.value))
        self._fig(self.fig_sv_prof, sv.plot_bidimensional_profile(res, figsize=(7, 1.7 * (len(vals) + 1) + 0.7)))
        self._say("Bi-dimensional profile: position along the filament against distance from its core.", "success")

    def _sv_compute_frames(self):
        from examm import streamviz as sv

        m, vals = self._sv_target()
        self._sv_frames = sv.comoving_frames(m, self.X, vals, side=float(self.w_sv_side.value), upsample=int(self.w_sv_up.value),
                                             half_thickness=float(self.w_sv_thick.value) or None)
        n = len(self._sv_frames["centers"])
        self.w_sv_frame.param.update(start=0, end=n - 1, value=n // 2, visible=True)
        self.dl_sv_movie.visible = True
        pts = self._sv_frames["n_points"]
        self.sv_info.object = (f"{n} frames; points on a plane: median {int(np.median(pts))}, fewest {int(pts.min())}."
                               + (" Few points per plane: increase the half thickness." if np.median(pts) < 20 else ""))
        self._sv_show_frame()
        self._say(f"Computed {n} co-moving frames; move the slider to walk along the filament.", "success")

    def _sv_show_frame(self):
        from examm import streamviz as sv

        if self._sv_frames is None:
            return
        self._fig(self.fig_sv_frame, sv.plot_comoving_frame(self._sv_frames, int(self.w_sv_frame.value), self.X))

    def _dl_frames_movie(self):
        from examm import streamviz as sv

        if self._sv_frames is None:
            raise ValueError("compute the co-moving frames first")
        with tempfile.TemporaryDirectory() as tmp:
            out = sv.animate_comoving_frames(self._sv_frames, Path(tmp) / "frames.mp4", self.X)
            self.dl_sv_movie.filename = "comoving_frames" + out.suffix
            return io.BytesIO(out.read_bytes())

    # ------------------------------------------------------------------ 3-D view
    def _view(self):
        feats = self.feature_names
        if self.w_pca.value and len(feats) > 2:
            return viz.View(self.X, projection="pca", names=feats)
        axes = [feats.index(w.value) for w in self.w_ax if w.value in feats]
        axes = list(dict.fromkeys(axes))
        if len(axes) < 2:
            axes = list(range(min(3, len(feats))))
        return viz.View(self.X, axes=axes, names=feats)

    def refresh_view(self):
        if self.X is None:
            return
        L, v = self.learner, self._view()
        on = set(self.w_layers.value)
        sc = Scene(labels=tuple(v.labels), flat=v.ndim == 2)
        rng = np.random.default_rng(0)

        def sub(n):
            return np.sort(rng.choice(n, self.max_display, replace=False)) if n > self.max_display else np.arange(n)

        focus = self.w_manifold.value if (L is not None and L._done >= 4 and self.w_manifold.value is not None) else -1
        if "Raw data" in on:
            s = sub(len(self.X))
            col = "#b9b8b2"
            if self.w_colour.value != "(none)":
                col = _colormap(self.table[:, self.columns.index(self.w_colour.value)])[s]
            sc.points("raw data", v(self.X[s]), col, 1.5, 0.6)
        if L is not None and L._done >= 1:
            pre = L.preprocess_
            s = sub(len(pre))
            if "Kept points (noisy)" in on:
                sc.points("kept points", v(pre.noisy[s]), "#8fb5e6", 1.5, 0.7)
            if "Diffused points" in on:
                sc.points("diffused points", v(pre.diffused[s]), "#eb6834", 1.8, 0.9)
            if "Dimensionality index" in on and L._done >= 2:
                idx = L.dim_result_.index[s]
                for j in np.unique(idx):
                    sc.points("undetermined" if j == 0 else f"{j}-D", v(pre.diffused[s][idx == j]),
                              viz.DIM_COLORS.get(int(j), "#0b0b0b"), 1.8)
        if L is not None and L._done >= 3:
            fitted = L._done >= 4
            items = list(enumerate(L.manifolds_)) if fitted else list(enumerate(g for gs in L.graphs_.values() for g in gs))
            for k, it in items:
                if focus >= 0 and k != focus:
                    continue
                g = it.graph if fitted else it
                c = viz.color_of(k)
                if "Attached points" in on:
                    pts = L.X_[it.point_indices] if fitted else L.preprocess_.noisy[g.points]
                    sc.points(f"points of {k}", v(pts[sub(len(pts))]), c, 1.6, 0.55)
                if "Graphs" in on:
                    sc.lines(f"graph {k}", v(g.nodes), g.edges, c, 2.5)
                if "AGTM (fitted graph)" in on and fitted and it.agtm is not None:
                    sc.lines(f"AGTM {k}", v(it.agtm.means_), g.edges, c if "Graphs" not in on else "#0b0b0b", 2.5)
        if "Scale r" in on and L is not None:
            ref = L.preprocess_.diffused if L._done >= 1 else self.X
            centre = ref[np.argmin(np.linalg.norm(ref - np.median(ref, 0), axis=1))]
            sc.sphere("scale r", v(centre[None])[0], float(L.r))
        if not sc.layers:
            sc.points("raw data", v(self.X[sub(len(self.X))]), "#b9b8b2", 1.5, 0.3)
        pane, used = render(sc, self.backend, pane=getattr(self, "_pane", None),
                            camera=(self.w_az.value, self.w_el.value, self.w_zoom.value))
        if pane is not getattr(self, "_pane", None):
            self._pane = pane
            self.view_holder[:] = [pane]
        shown = min(len(self.X), self.max_display)
        self.backend_note.object = (
            f"3-D renderer: {used}" + (" (drag to rotate, scroll to zoom)" if used == "plotly" else
                                       " (server-rendered image; use the sliders to move the camera)")
            + (f" \u00b7 showing {shown:,} of {len(self.X):,} points per layer" if len(self.X) > shown else ""))

    # ------------------------------------------------------------------ export
    def _need_fit(self):
        if self.learner is None or self.learner._done < 4:
            raise ValueError("run all steps first")
        return self.learner

    def _dl_model(self):
        return io.BytesIO(pickle.dumps(self._need_fit(), protocol=pickle.HIGHEST_PROTOCOL))

    def _dl_labels(self):
        L = self._need_fit()
        head = ",".join(["row", "kept", "dim_index", "manifold"] + self.feature_names)
        kept = np.zeros(len(self.X), int)
        kept[L.preprocess_.indices] = 1
        buf = io.StringIO()
        cols = np.column_stack([np.arange(len(self.X)), kept, L.dim_index_, L.labels_])
        np.savetxt(buf, np.column_stack([cols, self.table[:, [self.columns.index(f) for f in self.feature_names]]]),
                   delimiter=",", header=head, comments="", fmt=["%d"] * 4 + ["%.10g"] * len(self.feature_names))
        return io.BytesIO(buf.getvalue().encode())

    def _dl_viewer(self):
        """Standalone HTML viewer of the session: one global model built on the graphs that were found."""
        from examm.hierarchy import HierarchicalAGTM
        from examm.viewer import export_viewer

        L = self._need_fit()
        graphs = [g for gs in L.graphs_.values() for g in gs]
        if not graphs:
            raise ValueError("no graphs were found; nothing to show")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            H = HierarchicalAGTM(background=max(float(self.w_bg.value), 0.01),
                                 agtm_params=self._params()["agtm_params"]).set_data(self.X)
            H.add_children(graphs)
        H.detections_ = [{"level": 2, "parent": "root", "r": float(L.r), "learner": L, "region": np.arange(len(self.X))}]
        extra = {c: self.table[:, k] for k, c in enumerate(self.columns) if c not in self.feature_names}
        with tempfile.TemporaryDirectory() as tmp:
            out = export_viewer(Path(tmp) / "viewer.html", self.X, {"GUI session": H}, names=self.feature_names,
                                title="examm viewer", extra=extra)
            return io.BytesIO(out.read_bytes())

    def _dl_summary(self):
        import matplotlib.pyplot as plt

        v = self._view()
        fig = viz.plot_summary(self._need_fit(), names=self.feature_names,
                               **(dict(projection="pca") if v.basis is not None else dict(axes=v.axes)))
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=150)
        plt.close(fig)
        buf.seek(0)
        return buf

    def script(self) -> str:
        """Python script reproducing the current session."""
        p = self._params()
        feats = self.feature_names
        pre, idx = self._as_code(self._main_steps()), self._as_code(self._index_steps())
        if self.source in EXAMPLES:
            fn = {"Three interlocked toroids (3-D)": "make_three_toroids(scale=0.6, random_state=0)",
                  "Mixed manifolds D1 (3-D)": "make_d1(random_state=0)",
                  "Three filaments (3-D)": "make_filaments(random_state=0)",
                  "3-torus (4-D)": "make_torus3(random_state=0)"}[self.source]
            load = f"table = examm.datasets.{fn}.X\ncolumns = {self.columns!r}"
        else:
            load = f"table, columns = examm.io.load_table({self.source!r})  # adjust the path"
        lines = [
            "# Generated by examm-gui: reproduces the session.", "import numpy as np", "import examm",
            "from examm.preprocess import LAAT, MBMS, SAF, DensityFilter, DensityThin, Pipeline", "", load,
            f"features = {feats!r}", "X = table[:, [columns.index(f) for f in features]]",
        ]
        if self.w_standardise.value:
            lines.append("X = (X - X.mean(0)) / X.std(0)")
        lines += [
            "", f"preprocess = {pre}", f"index_preprocess = {idx}      # where the dimensionality index is estimated", "",
            "learner = examm.MultiManifoldLearner(",
            f"    r={p['r']!r}, preprocess=preprocess, smooth_factor={p['smooth_factor']!r},",
            f"    index_preprocess=index_preprocess, index_density_weight={p['index_density_weight']!r},",
            f"    graph={p['graph']!r}, dims={p['dims']!r}, eta={p['eta']!r}, beta={p['beta']!r}, "
            f"min_nodes={p['min_nodes']!r},",
            f"    epsilon={p['epsilon']!r}, background={p['background']!r},",
            f"    agtm_params={p['agtm_params']!r},".replace("inf", "np.inf"),
            "    random_state=0,", ").fit(X)", "", "print(learner.summary())",
            "# learner.manifolds_, learner.dim_index_, learner.labels_, learner.score_samples(X) ...",
            "# import examm.viz as viz; viz.plot_summary(learner).savefig('summary.png', dpi=150)", "",
        ]
        return "\n".join(lines)

    def _dl_script(self):
        if self.learner is None:
            raise ValueError("load a data set first")
        return io.BytesIO(self.script().encode())


def create_app(path: str | None = None, backend: str = "plotly", max_display: int = 40000):
    """Build the GUI and return the Panel object (servable; also usable inside a notebook)."""
    import panel as pn

    pn.extension("plotly", notifications=False)
    return App(path, backend, max_display).layout()


def main(argv=None):
    parser = argparse.ArgumentParser(prog="examm-gui", description="Step-by-step GUI for examm (opens in the browser).")
    parser.add_argument("file", nargs="?", help="table to load at start (csv, txt, npy, npz, h5, fits)")
    parser.add_argument("--port", type=int, default=5006)
    parser.add_argument("--address", default="localhost", help="address to listen on (default: localhost only)")
    parser.add_argument("--backend", choices=["plotly", "pyvista"], default="plotly",
                        help="3-D view: plotly = interactive WebGL in the browser (default); pyvista = "
                             "server-rendered images with camera sliders (needs off-screen OpenGL)")
    parser.add_argument("--max-display", type=int, default=40000, help="points drawn per layer")
    parser.add_argument("--no-browser", action="store_true", help="do not open a browser window")
    args = parser.parse_args(argv)
    try:
        import panel as pn
    except ImportError as exc:  # pragma: no cover
        raise SystemExit("The GUI needs extra packages: pip install 'examm[gui]'") from exc
    if args.backend == "pyvista" and not pyvista_available():
        raise SystemExit("PyVista cannot render off-screen on this machine; install an OpenGL library "
                         "(e.g. sudo apt install libosmesa6) or use the default --backend plotly")
    pn.extension("plotly")
    print(f"examm-gui: http://{args.address}:{args.port}  (Ctrl+C to stop)")
    pn.serve(lambda: create_app(args.file, args.backend, args.max_display), port=args.port, address=args.address,
             show=not args.no_browser, title="examm", websocket_origin="*" if args.address != "localhost" else None,
             threaded=False)


if __name__ == "__main__":
    main()
