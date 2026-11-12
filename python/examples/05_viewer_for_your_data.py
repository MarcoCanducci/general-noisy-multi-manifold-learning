"""Fit hierarchies on your own table and write the interactive viewer for them.

    python 05_viewer_for_your_data.py mydata.npy            # or .csv, .txt, .npz, .h5, .fits
    python 05_viewer_for_your_data.py mydata.csv --labels properties.csv -o my_viewer.html

Edit RUNS below to choose the scales and the preprocessing; every entry becomes one "version"
in the viewer. --labels adds per-point columns for colouring (they can also be loaded later in
the page itself).
"""
import argparse
import warnings

import numpy as np

import examm
from examm.viewer import export_viewer

parser = argparse.ArgumentParser()
parser.add_argument("data")
parser.add_argument("--labels", help="table with one row per point (same order) of extra labels")
parser.add_argument("-o", "--output", default="viewer.html")
args = parser.parse_args()

X, names = examm.io.load_table(args.data)

def plain():
    return examm.Pipeline([examm.SAF(), examm.DensityFilter(radius_factor=1.0, min_points=10)])

def thinned(radius, target=30):
    return examm.Pipeline([examm.DensityThin(radius=radius, target=target), examm.SAF(),
                           examm.DensityFilter(radius_factor=1.0, min_points=10)])

# name -> keyword arguments of examm.build_hierarchy (radii go from coarse to fine)
RUNS = {
    "Thinned hierarchy (r=0.7, then r=0.4)": dict(radii=[0.7, 0.4], preprocess=[thinned(0.4), plain()]),
    "One level, no thinning (r=0.4)": dict(radii=[0.4], preprocess=[plain()]),
}

versions = {}
with warnings.catch_warnings():
    warnings.simplefilter("ignore", examm.DegenerateModelWarning)
    for name, kw in RUNS.items():
        versions[name] = examm.build_hierarchy(X, background=0.01, agtm_params=dict(min_resp=3), random_state=0, **kw)
        print(name); print(versions[name].summary())

extra = None
if args.labels:
    import csv

    with open(args.labels, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if len(rows) != len(X):
        raise SystemExit(f"{args.labels} has {len(rows)} rows, the data {len(X)}")
    extra = {}
    for key in rows[0]:
        col = [row[key] for row in rows]
        try:
            extra[key] = np.array([float(v) if v != "" else np.nan for v in col])
        except ValueError:
            extra[key] = np.array(col)

out = export_viewer(args.output, X, versions, names=names, extra=extra, title=f"examm: {args.data}")
examm.io.save(versions, args.output.rsplit(".", 1)[0] + "_models.pkl")
print("wrote", out)
