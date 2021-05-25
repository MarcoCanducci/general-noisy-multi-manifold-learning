# general-noisy-multi-manifold-learning

This is the original matlab implementation of the paper:
"Probabilistic modelling of general noisy multi-manifolddata sets" by M. Canducci, P. Tino, M. Mastropietro, 2021

It is a complete methodology aimed at recovering low-dimensional noisy structures from noisy point clouds embedded in high dimensional uniform noise. The methodology is charachterized by different steps:

1. Diffusion via the SAF algorithm (check the corresponding paper for parameter tuning: https://www.cs.umd.edu/~zwicker/publications/StructureAwareDataConsolidation-PAMI2017.pdf)

2. Filtering based on difference of local dimensionalities between diffused and noisy data sets.

3. Dimensionality index estimation based on local PCA and clustering onto simplex of multimodal distributions. 

4. Multi-Manifold crawling and alternative graph construction for extraction of abstract and embedded graphs describing the geometry and topology of the detected manifolds

5. AGTM: A reformulation of Generative Topographic Mapping (https://www.microsoft.com/en-us/research/wp-content/uploads/1998/01/bishop-gtm-ncomp-98.pdf) 
that initializes the latent space as an abstract graphs and build manifold aligned noise models guaranteed to maximize the likelihood.

Additionally, we provide the codes for reproducing the data sets used in the paper for validation, with file Example.m that shows hot to apply the code to a data set.
