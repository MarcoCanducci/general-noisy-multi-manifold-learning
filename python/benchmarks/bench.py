"""Time every stage of the pipeline on the paper's data sets at increasing size.

    python benchmarks/bench.py            # default sizes
    python benchmarks/bench.py 1 4 10     # density multipliers
"""
import sys
import time
import warnings

import manydream
from manydream import datasets

scales = [float(a) for a in sys.argv[1:]] or [1.0, 3.0]
warnings.simplefilter("ignore")
print(f"threads: {manydream._core.get_num_threads()}, OpenMP: {manydream._core.has_openmp()}")
for name, make in (("three toroids", datasets.make_three_toroids), ("mixed manifolds D1", datasets.make_d1)):
    for s in scales:
        ds = make(scale=s, random_state=0)
        # the filter threshold scales with the sampling density
        pre = manydream.Pipeline([manydream.SAF(), manydream.DensityFilter(min_points=max(6, round(6 * s)))])
        L = manydream.MultiManifoldLearner(r=0.1, preprocess=pre, random_state=0)
        t = time.perf_counter()
        L.fit(ds.X)
        total = time.perf_counter() - t
        stages = "  ".join(f"{k} {v:6.2f}s" for k, v in L.timings_.items())
        print(f"{name:20s} n={len(ds.X):8,d}  manifolds={len(L.manifolds_):3d}  total {total:7.2f}s  |  {stages}")
