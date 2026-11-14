"""Three interlocked noisy toroids in uniform background noise (data set D2 of the paper)."""
import manydream
from manydream import datasets, viz

ds = datasets.make_three_toroids(random_state=0)
learner = manydream.MultiManifoldLearner(r=0.1, random_state=0, verbose=True).fit(ds.X)
print(learner.summary())

viz.plot_summary(learner).savefig("three_toroids_summary.png", dpi=150)
viz.plot_agtm(learner, ds.X).savefig("three_toroids_agtm.png", dpi=150)
