"""One-dimensional manifolds (filaments), the 1-DREAM use case."""
import examm
from examm import datasets, viz

ds = datasets.make_filaments(random_state=0)
learner = examm.MultiManifoldLearner(r=0.1, dims=[1], random_state=0).fit(ds.X)
print(learner.summary())
for row in viz.filament_table(learner):
    print(row)

viz.plot_filaments(learner, ds.X).savefig("filaments.png", dpi=150)
longest = max(learner.manifolds_, key=lambda m: m.graph.total_length())
viz.plot_profile(longest).savefig("filament_curvature.png", dpi=150)
viz.plot_profile(longest, ds.X, ds.X[:, 2], label="x3").savefig("filament_x3_profile.png", dpi=150)
viz.plot_cross_sections(longest, ds.X).savefig("filament_cross_sections.png", dpi=150)
