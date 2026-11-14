import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

import manydream  # noqa: E402
from manydream import datasets, viz  # noqa: E402


@pytest.fixture(scope="module")
def fitted():
    ds = datasets.make_filaments(random_state=0)
    extra = datasets.noisy_sample(datasets.torus(4000, random_state=1) + [2.2, 0.5, 0.5], 0.03, 6e4, 2)
    X = np.vstack([ds.X, extra])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        L = manydream.MultiManifoldLearner(r=0.1, random_state=0, agtm_params={"max_iter": 10}).fit(X)
    return X, L


def test_fixture_has_curves_and_a_surface(fitted):
    _, L = fitted
    dims = [m.dim for m in L.manifolds_]
    assert dims.count(1) >= 3 and dims.count(2) >= 1


def test_every_figure_renders(fitted, tmp_path):
    X, L = fitted
    surf = next(m for m in L.manifolds_ if m.dim == 2)
    fil = next(m for m in L.manifolds_ if m.dim == 1)
    figs = {
        "data": viz.plot_data(X), "data_c": viz.plot_data(X, color=X[:, 0], color_label="x"),
        "data_l": viz.plot_data(X, labels=np.maximum(L.labels_ + 1, 0)),
        "pre": viz.plot_preprocess(X, L.preprocess_), "counts": viz.plot_density_counts(X, 0.05, 6),
        "dim": viz.plot_dim_index(L.preprocess_.diffused, L.dim_result_), "simplex": viz.plot_simplex(L.dim_result_),
        "conf": viz.plot_confusion(L.dim_index_, (L.labels_ >= 0).astype(int)),
        "graphs": viz.plot_graphs(L, X, panels=True), "graphs2": viz.plot_graphs(L.graphs_),
        "agtm": viz.plot_agtm(L, X), "loglik": viz.plot_loglik(L), "curv": viz.plot_curvature(surf, X),
        "unfold": viz.plot_unfolded(surf, surf.agtm.edge_curvature(), label="curvature"),
        "nodes": viz.plot_node_values(surf, X, X[:, 2], label="z"), "fil": viz.plot_filaments(L, X),
        "prof": viz.plot_profile(fil), "prof2": viz.plot_profile(fil, X, X[:, 2], label="z"),
        "cross": viz.plot_cross_sections(fil, X), "summary": viz.plot_summary(L),
        "crawl": viz.plot_crawl_progress(datasets.circle(2000, 1, 2, 0), 0.2, 1, steps=(5, None)),
    }
    for name, fig in figs.items():
        out = tmp_path / f"{name}.png"
        fig.savefig(out, dpi=50)
        assert out.stat().st_size > 2000, name
    plt.close("all")
    rows = viz.filament_table(L)
    assert len(rows) >= 3 and all(r["length"] > 0 for r in rows)


def test_views_for_more_than_three_features(rng):
    X = datasets.torus3(3000, rng)
    for kw in (dict(axes=(0, 3)), dict(axes=(1, 2, 3)), dict(projection="pca"), dict(projection="pca", axes=(0, 1))):
        fig = viz.plot_data(X, names=["a", "b", "c", "e"], **kw)
        assert len(fig.axes) >= 1
    v = viz.View(X, projection="pca")
    assert v(X).shape == (3000, 3) and v.labels == ["PC1", "PC2", "PC3"]
    with pytest.raises(ValueError):
        viz.View(X, axes=(0, 9))
    plt.close("all")


def test_chain_order_follows_the_filament():
    t = np.linspace(0, 1, 12)
    perm = np.random.default_rng(0).permutation(12)
    nodes = np.c_[t, t**2][perm]
    inv = np.argsort(perm)
    g = manydream.ManifoldGraph(nodes, [[inv[i], inv[i + 1]] for i in range(11)], 1)
    order, s = viz.chain_order(g)
    assert np.allclose(np.abs(np.diff(nodes[order][:, 0])), t[1] - t[0]) and np.all(np.diff(s) > 0)
