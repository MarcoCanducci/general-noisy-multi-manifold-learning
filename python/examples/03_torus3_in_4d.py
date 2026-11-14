"""A 3-dimensional manifold (3-torus) embedded in 4-D (paper section 7)."""
import manydream
from manydream import datasets, viz

X = datasets.torus3(grid=30)                      # regular sampling of the 3-torus, 27 000 points
graph = manydream.crawl(X, r=1.0, dim=3, seed=0)
print(graph, "| maximum degree", graph.degree.max())
viz.plot_graphs([graph], X, axes=(0, 1, 2)).savefig("torus3_dims123.png", dpi=150)
viz.plot_graphs([graph], X, axes=(1, 2, 3)).savefig("torus3_dims234.png", dpi=150)
