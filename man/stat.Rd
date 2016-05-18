\name{stat}
\alias{stat}
\title{Fourthcorner analysis for transcriptomic data}
\description{
The fourthcorner analysis tests for significant associations between each sample covariate and each gene covariate by statistical permutation tests. The sample and gene covariates can be categorical and/or quantitative.

The input has to be given as dataframe or matrix. Dataframe/matrix L [n x p] contains transcriptomic data of p samples across n genes, dataframe/matrix R [n x m] contains m gene covariates across the n genes and dataframe/matrix Q [p x s] contains s sample covariates across the p samples. Alternatively, objects of the class ExpressionSet (with assayData, phenoData and featureData) can be used as input. If the argument ExprSet is missing, the function will use the dataframes/matrices R, L and Q as input.

The number of permutations is set to 9999 per default to assure significance of p-values after multiple testing correction. As computation time increases with size of the matrices/dataframes and with number of permutations, parallelization across multiple cores is highly recommended. Per default, all except one CPU cores on the current host are used.

Genes can be filtered with respect to their expression variance before analysis (argument exprvar); the function will automatically discard the gene covariates which do not annotate any of the remaining genes.

Warning: If R and Q are given as matrices, they will be converted to dataframes at the beginning of the function.

Warning: If R or Q is missing, it will be replaced by an identity matrix.
}
\usage{
stat(ExprSet, R=NULL, L=NULL, Q=NULL, npermut=9999, padjust="BH",
         nrcor=detectCores()-1, exprvar=1)
}

\arguments{
  \item{ExprSet}{An ExpressionSet of the \emph{Biobase} package.
  The ExpressionSet is used as default input. If no ExpressionSet is given, the individual dataframes/matrices R, L and Q can be used as input.}
  \item{R}{A dataframe/matrix containing information about each gene. The number of rows in R must match the number of rows in L. If R is missing, it will be replaced by an identity matrix [n x n].}
  \item{L}{A dataframe/matrix of gene expression values of genes across samples.}
  \item{Q}{A dataframe/matrix containing information about each sample. The number of rows in Q must match the number of columns in L. If Q is missing, it will be replaced by an identity matrix [p x p].}
  \item{npermut}{The number of permutations.}
  \item{padjust}{The method of multiple testing adjustment of the pvalues, see p.adjust.methods for all methods implemented in R.}
  \item{nrcor}{The number of cores to be used.}
  \item{exprvar}{The fraction of most variably expressed genes to take into account. If the functions `stat` and `ord` shall be combined, this value has to be the same in both analyses.}
}

\details{
Dependent on the covariate combination, a statistic is calculated based on matrix multiplication of the three tables. This statistic amounts to a correlation coefficient for the association between quantitative-quantitative and quantitative-categorical variables and to a Chi2-related statistic for the association between categorical-categorical variables.
}
\value{
The function returns a list of class stat where:
  \item{stat}{is a cross table (m x s) with the values of the original statistical tests per covariate combination.}
  \item{pvalue, adj.pvalue}{are cross tables (m x s) which contain the p-values and adjusted p-values, respectively, of the permutation tests per covariate combination.}
  \item{adjust.method}{shows the applied multiple testing adjustment method.}
  \item{npermut}{gives the number of permutations per permutation test.}
  \item{ngenes}{gives the number of analysed genes ("all" in the case of no filtering of the genes).}
  \item{call}{gives the original call of the function.}
}
\author{Lara Urban}
\examples{
data(Baca)
statBaca <- stat(ExprSet = Baca, npermut = 999, padjust = "BH", nrcor = 2, exprvar = 1)
statBaca$adj.pvalue
plot(statBaca)
}
