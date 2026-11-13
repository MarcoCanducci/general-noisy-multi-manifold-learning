"""Manifold crawling (paper sections 3.2 and 4.2) and the graph container shared by all stages."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import sparse
from scipy.sparse import csgraph

from ._util import as_points, build_tree, check_radius, rng_from

__all__ = ["ManifoldGraph", "crawl", "multi_manifold_crawl", "attach_points"]


@dataclass(eq=False)
class ManifoldGraph:
    """Skeleton of one manifold: an abstract graph and its embedding in data space.

    Attributes
    ----------
    nodes : (K, d) array
        Embedded node positions.
    edges : (E, 2) int array
        Undirected edges between node indices.
    dim : int
        Intrinsic dimension j of the manifold.
    frames : (K, j, d) array or None
        Orthonormal tangent frame at every node (rows span the local tangent space).
    node_index : (K,) int array or None
        Row of each node in the point set that was crawled.
    points : (n,) int array or None
        Rows of the noisy data set attached to this manifold.
    """

    nodes: np.ndarray
    edges: np.ndarray
    dim: int
    frames: np.ndarray | None = None
    node_index: np.ndarray | None = None
    points: np.ndarray | None = None
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        self.nodes = np.ascontiguousarray(self.nodes, dtype=float)
        self.edges = np.asarray(self.edges, dtype=np.int64).reshape(-1, 2)

    @property
    def n_nodes(self) -> int:
        return len(self.nodes)

    @property
    def n_edges(self) -> int:
        return len(self.edges)

    @property
    def edge_lengths(self) -> np.ndarray:
        return np.linalg.norm(self.nodes[self.edges[:, 0]] - self.nodes[self.edges[:, 1]], axis=1)

    @property
    def degree(self) -> np.ndarray:
        return np.bincount(self.edges.ravel(), minlength=self.n_nodes)

    def adjacency(self, weighted: bool = False) -> sparse.csr_matrix:
        """Symmetric sparse adjacency matrix (edge lengths as weights if ``weighted``)."""
        w = self.edge_lengths if weighted else np.ones(self.n_edges)
        i, j = self.edges.T
        A = sparse.coo_matrix((np.r_[w, w], (np.r_[i, j], np.r_[j, i])), shape=(self.n_nodes,) * 2)
        return A.tocsr()

    def components(self) -> list["ManifoldGraph"]:
        """Split into connected components."""
        n, lab = csgraph.connected_components(self.adjacency(), directed=False)
        if n <= 1:
            return [self]
        out = []
        for c in range(n):
            keep = np.flatnonzero(lab == c)
            remap = -np.ones(self.n_nodes, dtype=np.int64)
            remap[keep] = np.arange(len(keep))
            e = self.edges[(lab[self.edges[:, 0]] == c)]
            out.append(ManifoldGraph(self.nodes[keep], remap[e], self.dim,
                                     None if self.frames is None else self.frames[keep],
                                     None if self.node_index is None else self.node_index[keep],
                                     meta=dict(self.meta)))
        return out

    def total_length(self) -> float:
        return float(self.edge_lengths.sum())

    def to_networkx(self):
        """Return a ``networkx.Graph`` with node attribute ``pos`` and edge attribute ``length``."""
        import networkx as nx

        G = nx.Graph()
        for k, p in enumerate(self.nodes):
            G.add_node(k, pos=p)
        for (a, b), w in zip(self.edges, self.edge_lengths):
            G.add_edge(int(a), int(b), length=float(w))
        return G

    def __repr__(self) -> str:
        n = "?" if self.points is None else len(self.points)
        return f"ManifoldGraph(dim={self.dim}, nodes={self.n_nodes}, edges={self.n_edges}, points={n})"


def crawl(X, r: float, dim: int, *, seed: int | None = None, eta: float = 0.75, beta: float = 0.4,
          min_points: int = 5, max_nodes: int | None = None, random_state=None, tree=None,
          alive=None) -> ManifoldGraph:
    """Crawl a single manifold of intrinsic dimension ``dim`` through the point cloud ``X``.

    Starting from ``seed`` the graph grows along the local tangent spaces: each node proposes
    ``2 * dim`` candidates at distance ``eta * r`` along its (transported) tangent frame; a
    candidate within ``beta * r`` of an existing node is merged with it, otherwise it becomes a
    new node. Every node has at most ``2 * dim`` edges.

    ``alive`` (optional uint8 mask, modified in place) restricts the points that may become nodes.
    """
    X = as_points(X)
    r = check_radius(r)
    n, d = X.shape
    if not 1 <= dim <= d:
        raise ValueError(f"dim must be between 1 and {d}")
    if not (0 < beta < eta <= 1):
        raise ValueError("parameters must satisfy 0 < beta < eta <= 1")
    tree = tree if tree is not None else build_tree(X)
    if alive is None:
        alive = np.ones(n, dtype=np.uint8)
    if seed is None:
        cand = np.flatnonzero(alive)
        if len(cand) == 0:
            raise ValueError("no points left to crawl")
        seed = int(rng_from(random_state).choice(cand))
    node_of = np.full(n, -1, dtype=np.int64)
    idx, edges, frames = tree.crawl(int(seed), alive, node_of, r, int(dim), float(eta), float(beta),
                                    int(min_points), int(max_nodes or np.iinfo(np.int64).max))
    return ManifoldGraph(X[idx], edges, int(dim), frames, idx)


def attach_points(graph: ManifoldGraph, noisy, r: float, tree=None, restrict=None) -> np.ndarray:
    """Rows of ``noisy`` lying within ``r`` of a node of ``graph`` (optionally among ``restrict``)."""
    tree = tree if tree is not None else build_tree(as_points(noisy))
    mask = tree.radius_union(graph.nodes, r).astype(bool)
    if restrict is not None:
        mask &= restrict
    return np.flatnonzero(mask)


def multi_manifold_crawl(diffused, r: float, dim: int, *, noisy=None, eta: float = 0.75,
                         beta: float = 0.4, min_points: int = 5, min_nodes: int | None = None,
                         min_residual: int = 20, random_state=None) -> list[ManifoldGraph]:
    """Extract all manifolds of intrinsic dimension ``dim`` from a diffused point set.

    Crawling is restarted from random seeds on the points not yet explained until fewer than
    ``min_residual`` remain. After each crawl the points within ``r`` of the new graph are
    removed. Graphs with fewer than ``min_nodes`` (default ``4 * dim + 4``) nodes are discarded.
    If ``noisy`` (row-aligned with ``diffused``) is given, each graph gets the indices of the
    noisy points within ``r`` of its nodes in ``graph.points``.
    """
    D = as_points(diffused)
    r = check_radius(r)
    n = len(D)
    if min_nodes is None:
        min_nodes = 4 * dim + 4
    rng = rng_from(random_state)
    tree = build_tree(D)
    alive = np.ones(n, dtype=np.uint8)
    node_of = np.full(n, -1, dtype=np.int64)
    graphs: list[ManifoldGraph] = []
    remaining = n
    big = np.iinfo(np.int64).max
    for seed in rng.permutation(n):
        if remaining < min_residual:
            break
        if not alive[seed]:
            continue
        idx, edges, frames = tree.crawl(int(seed), alive, node_of, r, int(dim), float(eta), float(beta),
                                        int(min_points), big)
        if len(idx) >= min_nodes:
            alive[tree.radius_union(D[idx], r).astype(bool)] = 0
            graphs.append(ManifoldGraph(D[idx], edges, int(dim), frames, idx))
        remaining = int(alive.sum()) if len(idx) else remaining - 1
    graphs.sort(key=lambda g: -g.n_nodes)
    if noisy is not None:
        Q = as_points(noisy, "noisy")
        if len(Q) != n:
            raise ValueError("noisy and diffused must be row-aligned")
        qt = build_tree(Q)
        for g in graphs:
            g.points = attach_points(g, Q, r, qt)
    return graphs
