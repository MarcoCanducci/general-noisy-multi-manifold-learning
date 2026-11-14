"""Results in the form they are used: skeletons, membership probabilities, validation.

* :func:`skeleton` / :func:`save_skeleton` – the nodes and edges of every fitted structure;
* :func:`membership` / :func:`save_membership` – for every point, the probability of belonging
  to each structure or to none of them (the background);
* :func:`validate` – comparison with known labels: how well noise is told from structure
  (false positives and negatives), and how well each true manifold is recovered.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ._util import as_points
from .agtm import MultiManifoldModel

__all__ = ["skeleton", "save_skeleton", "membership", "save_membership", "validate", "ValidationReport"]


def _fitted(learner):
    ms = [m for m in getattr(learner, "manifolds_", []) if m.agtm is not None]
    if not ms:
        raise ValueError("the learner has no fitted manifolds; run it to the 'agtm' stage first")
    return ms


def skeleton(learner, fitted: bool = True) -> dict:
    """Nodes and edges of every structure.

    ``fitted=True`` (default) gives the node positions after the model was trained, which lie
    on the mean curve or surface of the noisy data; ``False`` gives the crawled graph.
    Returns ``nodes`` (structured array: structure, node, dim, then one field per coordinate)
    and ``edges`` (structure, node_a, node_b, length).
    """
    ms = _fitted(learner)
    d = ms[0].graph.nodes.shape[1]
    names = [f"x{k + 1}" for k in range(d)]
    nrows, erows = [], []
    for k, m in enumerate(ms):
        P = m.agtm.means_ if fitted else m.graph.nodes
        E = np.asarray(m.graph.edges).reshape(-1, 2)
        nrows += [(k, i, m.dim, *p) for i, p in enumerate(P)]
        erows += [(k, int(a), int(b), float(np.linalg.norm(P[a] - P[b]))) for a, b in E]
    nodes = np.array(nrows, dtype=[("structure", int), ("node", int), ("dim", int)] + [(n, float) for n in names])
    edges = np.array(erows, dtype=[("structure", int), ("node_a", int), ("node_b", int), ("length", float)])
    return {"nodes": nodes, "edges": edges}


def _write(path, arr):
    cols = arr.dtype.names
    fmt = ["%d" if arr.dtype[c].kind in "iu" else "%.8g" for c in cols]
    np.savetxt(path, np.column_stack([arr[c] for c in cols]), delimiter=",", header=",".join(cols), comments="", fmt=fmt)


def save_skeleton(learner, prefix, fitted: bool = True) -> tuple:
    """Write ``<prefix>_nodes.csv`` and ``<prefix>_edges.csv``; returns the two paths."""
    sk = skeleton(learner, fitted)
    a, b = Path(f"{prefix}_nodes.csv"), Path(f"{prefix}_edges.csv")
    _write(a, sk["nodes"])
    _write(b, sk["edges"])
    return a, b


def membership(learner, X=None, background: float | None = None) -> dict:
    """Probability that each point belongs to each structure, or to none.

    The structures form a mixture, weighted by their numbers of points, with a uniform
    background over the extent of the data. ``background`` is its prior weight; by default the
    learner's own value, or 0.05 if that is zero (without a background every point is forced
    onto some structure). Points outside the extent of the training data are background. Returns ``proba`` (n, structures + 1; the last column is the
    background), ``structure`` (most probable one, -1 for background), ``on_structure``
    (1 - background probability) and ``log_density``.
    """
    ms = _fitted(learner)
    X = learner.X_ if X is None else as_points(X)
    bg = float(learner.background if background is None else background) or 0.05
    lo, hi = learner.X_.min(0), learner.X_.max(0)
    model = MultiManifoldModel([m.agtm for m in ms], weights=[m.n_points for m in ms], background=bg, bounds=(lo, hi))
    P = model.predict_proba(X)
    P = np.where(np.isfinite(P), P, 0.0)
    outside = ~((X >= lo) & (X <= hi)).all(1)          # beyond the extent of the data: no structure can claim the point
    P[outside] = 0.0
    P[outside, -1] = 1.0
    lab = np.argmax(P, axis=1)
    lab[lab == len(ms)] = -1
    return {"proba": P, "structure": lab, "on_structure": 1.0 - P[:, -1], "log_density": model.score_samples(X),
            "background": bg}


def save_membership(learner, path, X=None, background: float | None = None) -> Path:
    """Write the membership table as CSV: one line per point."""
    r = membership(learner, X, background)
    K = r["proba"].shape[1] - 1
    head = ["row", "structure", "p_on_structure", "p_background"] + [f"p_structure_{k}" for k in range(K)]
    T = np.column_stack([np.arange(len(r["structure"])), r["structure"], r["on_structure"], r["proba"][:, -1], r["proba"][:, :-1]])
    np.savetxt(path, T, delimiter=",", header=",".join(head), comments="", fmt=["%d", "%d"] + ["%.6g"] * (K + 2))
    return Path(path)


@dataclass
class ValidationReport:
    """Outcome of :func:`validate`."""

    detection: dict                       # noise against structure: counts, precision, recall
    manifolds: list = field(default_factory=list)   # one dict per true manifold
    spurious: list = field(default_factory=list)    # predicted structures made mostly of noise

    def __str__(self) -> str:
        d = self.detection
        out = ["Noise against structure",
               f"  true positives  {d['tp']:>8,d}   points of a manifold assigned to a structure",
               f"  false negatives {d['fn']:>8,d}   points of a manifold left out",
               f"  false positives {d['fp']:>8,d}   noise points assigned to a structure",
               f"  true negatives  {d['tn']:>8,d}   noise points left out",
               f"  precision {d['precision']:.3f}   recall {d['recall']:.3f}   F1 {d['f1']:.3f}", "",
               "Each true manifold                 points   recall  purity  pieces  matched structure"]
        for m in self.manifolds:
            out.append(f"  {str(m['name'])[:28]:28s} {m['points']:>8,d}   {m['recall']:6.3f}  {m['purity']:6.3f}  {m['pieces']:>6d}  "
                       + ("none" if m["structure"] < 0 else str(m["structure"])))
        if self.spurious:
            out += ["", "Structures made mostly of noise: " + ", ".join(f"{s['structure']} ({s['points']:,d} points, {s['noise_share']:.0%} noise)"
                                                                      for s in self.spurious)]
        return "\n".join(out)


def validate(predicted, truth, names=None, min_piece: float = 0.05) -> ValidationReport:
    """Compare an assignment of points to structures with known labels.

    Parameters
    ----------
    predicted : (n,) int
        Structure of every point, -1 (or any negative value) for points left out, e.g.
        ``learner.labels_`` or ``membership(learner)["structure"]``.
    truth : (n,) int
        True label of every point: 0 for noise, ``k >= 1`` for the k-th manifold (as in
        :class:`manydream.datasets.Dataset`).
    names : list, optional
        Names of the true labels (index 0 is the noise).
    min_piece : float
        A predicted structure counts as a piece of a true manifold if it holds at least this
        share of the manifold's points.

    The report gives, for noise against structure, the four counts with precision and recall;
    and for every true manifold its best-matching structure, the share of the manifold found
    there (recall), the share of that structure made of the manifold (purity), and into how
    many pieces the manifold was split.
    """
    pred, truth = np.asarray(predicted).astype(int), np.asarray(truth).astype(int)
    if pred.shape != truth.shape:
        raise ValueError("predicted and truth must have the same length")
    on, real = pred >= 0, truth > 0
    tp, fp, fn, tn = int((on & real).sum()), int((on & ~real).sum()), int((~on & real).sum()), int((~on & ~real).sum())
    prec = tp / max(tp + fp, 1)
    rec = tp / max(tp + fn, 1)
    det = dict(tp=tp, fp=fp, fn=fn, tn=tn, precision=prec, recall=rec, f1=2 * prec * rec / max(prec + rec, 1e-300))
    structures = np.unique(pred[on])
    size = {int(s): int((pred == s).sum()) for s in structures}
    rows = []
    for k in np.unique(truth[real]):
        mine = truth == k
        n = int(mine.sum())
        share = {int(s): int((mine & (pred == s)).sum()) for s in structures}
        best = max(share, key=share.get) if share and max(share.values()) > 0 else -1
        rows.append(dict(label=int(k), name=names[int(k)] if names is not None else f"manifold {int(k)}", points=n, structure=best,
                         recall=share[best] / n if best >= 0 else 0.0, purity=share[best] / size[best] if best >= 0 else 0.0,
                         pieces=sum(1 for v in share.values() if v >= min_piece * n), found=sum(share.values()) / n))
    spurious = [dict(structure=s, points=size[s], noise_share=float(((pred == s) & ~real).sum() / size[s]))
                for s in size if ((pred == s) & ~real).sum() > 0.5 * size[s]]
    return ValidationReport(det, rows, spurious)
