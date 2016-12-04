# covRNA

The R package covRNA provides a convenient and fast interface for testing and visualizing the relationship between sample and gene covariates mediated by gene expression data. The fourthcorner analysis tests the statistical significance of sample and gene covariate relationships by permutation tests while the RLQ visualizes relationships within and between sample and gene covariates.

The method is based on the powerful fourthcorner and RLQ analyses used in ecological research for the analysis of species abundance data. We have modified the algorithms of the R package ade4 to make the method suitable to the distributional characteristics of RNA-Seq read counts and microarray intensities and to provide a parallelized high-performance implementation to allow for the analysis of large-scale transcriptomics data. We further supply the user with methods for unsupervised gene filtering and plotting functions to ease the analysis workflow. cov-RNA can be easily applied to microarray or RNA-Seq data of any organism. We provide a vignette with an exemplified analysis of an RNA-Seq dataset of Bacillus anthracis.

To cite covRNA, please use citation("covRNA"). For further details, please refer to the vignette and the man pages.


