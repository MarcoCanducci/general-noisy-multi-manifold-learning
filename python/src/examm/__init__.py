"""examm: probabilistic modelling of general noisy multi-manifold data sets.

Python/C++ implementation of M. Canducci, P. Tino, M. Mastropietro,
"Probabilistic modelling of general noisy multi-manifold data sets",
Artificial Intelligence 302 (2022) 103579.
"""
from . import datasets, preprocess
from ._util import set_num_threads
from .dimension import DimIndexResult, dimensionality_index
from .preprocess import SAF, DensityFilter, Pipeline, Preprocessor, PreprocessResult

__version__ = "0.1.0"
