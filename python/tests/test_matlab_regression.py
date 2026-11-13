"""Soft regression against the workspace saved by the original MATLAB pipeline.

Exact equality is not expected (the Python package is a redesign); the checks are agreement rates
with the stored MATLAB results on the same input.
"""
import warnings
from pathlib import Path

import numpy as np
import pytest

import manydream

MAT = Path(__file__).resolve().parents[2] / "3ToroidsNoisy_BestModel.mat"


@pytest.mark.skipif(not MAT.is_file(), reason="MATLAB workspace not available")
def test_agreement_with_matlab_results_on_three_toroids():
    from scipy.io import loadmat

    m = loadmat(MAT, variable_names=["DNN", "Labels", "Idx"])
    X, kept_m, idx_m = np.ascontiguousarray(m["DNN"][:, :3]), m["Labels"].ravel() == 1, m["Idx"].ravel()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        L = manydream.MultiManifoldLearner(r=0.2, random_state=0, agtm_params={"max_iter": 10}).fit(X)
    kept = np.zeros(len(X), bool)
    kept[L.preprocess_.indices] = True
    assert (kept == kept_m).mean() > 0.95                       # density filter
    full = np.zeros(len(X), int)
    full[kept_m] = idx_m
    both = kept & kept_m
    assert (L.dim_index_[both] == full[both]).mean() > 0.97     # dimensionality index
    # MATLAB assigned 5083 / 5256 / 5075 points to the three toroids
    assert [mf.dim for mf in L.manifolds_] == [2, 2, 2]
    assert all(4500 < mf.n_points < 5600 for mf in L.manifolds_)
