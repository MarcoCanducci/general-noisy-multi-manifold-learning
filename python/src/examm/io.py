"""Reading tabular data and saving / loading fitted results."""
from __future__ import annotations

import pickle
from pathlib import Path

import numpy as np

__all__ = ["load_table", "save", "load"]


def load_table(path, key: str | None = None):
    """Load a table of numbers. Returns ``(array (n, c), column_names)``.

    Supported: ``.csv`` / ``.txt`` / ``.dat`` / ``.tsv`` (with or without a header line),
    ``.npy``, ``.npz``, ``.h5`` / ``.hdf5`` (needs h5py) and ``.fits`` tables (needs astropy).
    ``key`` selects the array / dataset / HDU where a file holds several.
    """
    path = Path(path)
    ext = path.suffix.lower()
    names = None
    if ext == ".npy":
        A = np.load(path)
    elif ext == ".npz":
        with np.load(path) as z:
            A = z[key or z.files[0]]
    elif ext in (".h5", ".hdf5"):
        import h5py

        with h5py.File(path, "r") as f:
            if key is None:
                found = []
                f.visititems(lambda n, o: found.append(n) if isinstance(o, h5py.Dataset) and o.ndim == 2 else None)
                if not found:
                    raise ValueError("no 2-D dataset found in the HDF5 file; pass key=")
                key = found[0]
            A = f[key][...]
    elif ext in (".fits", ".fit"):
        from astropy.table import Table

        t = Table.read(path, hdu=key) if key is not None else Table.read(path)
        names = [c for c in t.colnames if np.issubdtype(t[c].dtype, np.number) and t[c].ndim == 1]
        A = np.column_stack([np.asarray(t[c], float) for c in names])
    else:
        with open(path) as fh:
            first = fh.readline()
        delim = "," if "," in first else ("\t" if "\t" in first else None)
        tokens = [t.strip().strip('"') for t in (first.split(delim) if delim else first.split())]
        header = False
        try:
            [float(t) for t in tokens]
        except ValueError:
            header = True
        A = np.loadtxt(path, delimiter=delim, skiprows=1 if header else 0, ndmin=2, comments="#")
        if header:
            names = [t.lstrip("#").strip() for t in tokens][: A.shape[1]]
    A = np.asarray(A, dtype=float)
    if A.ndim == 1:
        A = A[:, None]
    if A.ndim != 2:
        raise ValueError(f"expected a 2-D table, got an array of shape {A.shape}")
    if names is None or len(names) != A.shape[1]:
        names = [f"x{k + 1}" for k in range(A.shape[1])]
    return A, names


def save(obj, path) -> None:
    """Save a fitted object (learner, AGTM, graph, ...) to ``path`` with pickle."""
    with open(path, "wb") as fh:
        pickle.dump(obj, fh, protocol=pickle.HIGHEST_PROTOCOL)


def load(path):
    """Load an object stored with :func:`save`. Only open files you trust."""
    with open(path, "rb") as fh:
        return pickle.load(fh)
