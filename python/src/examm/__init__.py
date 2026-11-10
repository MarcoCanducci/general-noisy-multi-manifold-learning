"""examm: probabilistic modelling of general noisy multi-manifold data sets.

Python/C++ implementation of M. Canducci, P. Tino, M. Mastropietro,
"Probabilistic modelling of general noisy multi-manifold data sets",
Artificial Intelligence 302 (2022) 103579.

Quick start::

    import examm
    learner = examm.MultiManifoldLearner(r=0.1, random_state=0).fit(X)
    print(learner.summary())
"""
from . import datasets, io, preprocess
from ._util import set_num_threads
from .agtm import AGTM, DegenerateModelError, DegenerateModelWarning, MultiManifoldModel
from .crawling import ManifoldGraph, crawl, multi_manifold_crawl
from .dimension import DimIndexResult, dimensionality_index
from .epsgraph import epsilon_graph
from .pipeline import STAGES, Manifold, MultiManifoldLearner
from .preprocess import SAF, DensityFilter, Pipeline, Preprocessor, PreprocessResult

__version__ = "0.1.0"

__all__ = [
    "AGTM", "DegenerateModelError", "DegenerateModelWarning", "DensityFilter", "DimIndexResult",
    "Manifold", "ManifoldGraph", "MultiManifoldLearner", "MultiManifoldModel", "Pipeline",
    "PreprocessResult", "Preprocessor", "SAF", "STAGES", "crawl", "datasets", "dimensionality_index",
    "epsilon_graph", "io", "multi_manifold_crawl", "preprocess", "set_num_threads",
]
