# covRNA

covRNA (covariate analysis of RNA-Seq data) is a fast and user-friendly R package which implements fourthcorner analysis and RLQ of transcriptomic data. 

Gene expression data normally comes with covariates of the samples and of the genes. To analyze associations between sample and gene covariates, the fourthcorner analysis tests the statistical significance of the associations by permutation tests while the RLQ visualizes associations within and be-tween the covariates.

The fourthcorner analysis and RLQ implemented in the ade4 package are adapted to easily analyze large-scale transcriptomic data. (1) Runtime and storage space are significantly reduced, (2) the analysis accounts for tran-scriptome-specific shapes of the empirical permutation distributions, (3) the analysis is rendered user-friendly by supplying automation, simple design-ing of plots and unsupervised gene filtering.

To cite covRNA, please use citation("covRNA"). For further details, please refer to the vignette and the man pages.
