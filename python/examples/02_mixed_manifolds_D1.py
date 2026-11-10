"""Curves, surfaces and a blob in heavy background noise (data set D1 of the paper), step by step."""
import numpy as np

import examm
from examm import datasets, viz

ds = datasets.make_d1(random_state=0)
X, r = ds.X, 0.1

# 1. denoise: SAF diffusion + density filter (any examm.preprocess.Preprocessor can be used instead)
pre = examm.Pipeline([examm.SAF(), examm.DensityFilter()]).fit_transform(X, r)
print(f"kept {len(pre)} of {len(X)} points")

# 2. dimensionality index of the diffused points
dim = examm.dimensionality_index(pre.diffused, r)
print("points per dimension:", dim.counts_per_dim())
viz.plot_confusion(dim.index, ds.dim[pre.indices]).savefig("d1_confusion.png", dpi=150)

# 3. crawl the 2-D manifolds and 4. fit an Abstract GTM on each
sel = dim.index == 2
graphs = examm.multi_manifold_crawl(pre.diffused[sel], r, dim=2, noisy=pre.noisy[sel], random_state=0)
models = [examm.AGTM(g).fit(pre.noisy[sel][g.points]) for g in graphs]
for g, m in zip(graphs, models):
    print(g, "->", m, "| median curvature", round(float(np.median(m.edge_curvature())), 2))

# 5. one density model for everything
mixture = examm.MultiManifoldModel(models)
print("mean log-likelihood of the 2-D points:", round(mixture.score(pre.noisy[sel]), 3))
viz.plot_graphs(graphs, X, panels=True).savefig("d1_graphs.png", dpi=150)
