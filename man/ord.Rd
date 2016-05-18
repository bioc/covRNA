\name{ord}
\alias{ord}
\title{RLQ for transcriptomic data}
\description{
The RLQ visualises the association between and within sample and gene covariates by ordination. It applies generalized singular value decomposition (GSVD) to the fourthcorner matrix, which contains the associations between the sample and gene covariates. This is realised by eigendecomposition of the covariance matrices of the fourthcorner matrix. The name RLQ refers to the three dataframes R, L and Q to be analyzed. 
The function `ord` automates the `rlq` function of the `ade4` package.

The input has to be given as dataframe or matrix. Dataframe/matrix L [n x p] contains transcriptomic data of p samples across n genes, dataframe/matrix R [n x m] contains m gene covariates across the n genes and dataframe/matrix Q [p x s] contains s sample covariates across the p samples. Alternatively, objects of the class ExpressionSet (with assayData, phenoData and featureData) can be used as input. If the argument ExprSet is missing, the function will use the dataframes/matrices R, L and Q as input.

Genes can be filtered with respect to their expression variance before analysis (argument exprvar); the function will automatically discard the gene covariates which do not annotate any of the remaining genes.

Warning: If R and Q are given as matrices, they will be converted to dataframes at the beginning of the function.

Warning: If R or Q is missing, it will be replaced by an identity matrix. Then, a principal component analysis of this matrix will be performed what might be time-consuming, depending on the size of the identity matrix.
}
\usage{
ord(ExprSet, R=NULL, L=NULL, Q=NULL, exprvar=1, nf=2)
}
\arguments{
  \item{ExprSet}{An ExpressionSet of the \emph{Biobase} package.
The ExpressionSet is used as default input. If no ExpressionSet is given, the individual dataframes/matrices R, L and Q can be used as input.}
  \item{R}{A dataframe/matrix containing information about each gene. The number of rows in R must match the number of rows in L. If R is missing, it will be replaced by an identity matrix [n x n].}
  \item{L}{A dataframe/matrix of gene expression values of genes across samples.}
  \item{Q}{A dataframe/matrix containing information about each sample. The number of rows in Q must match the number of columns in L. If Q is missing, it will be replaced by an identity matrix [p x p].}
  \item{exprvar}{The fraction of most variably expressed genes to take into account. If the functions `stat` and `ord` shall be combined, this value has to be the same in both analyses.}
  \item{nf}{The number of axes to be considered by ordination.}
}
\details{
The function automates the following steps. Firstly, Correspondence Analysis is applied to gene expression table L. Either Principal Component Analysis (only quantitative variables), Multiple Correspondence Analysis (only categorical variables) or Hillsmith analysis (quantitative and categorical variables) are applied to the covariate tables R and Q. Secondly, RLQ is applied to the results of these ordination methods.
}
\value{
The function returns a list ob class ord where:
  \item{call}{gives the original call of the function.}
  \item{rank}{gives the rank.}
  \item{nf}{gives number of axes to be considered by ordination.}
  \item{RV}{gives the RV coefficient.}
  \item{eig}{gives a vector of the eigenvalues.}
  \item{variance}{gives the variance explained by the axes.}
  \item{lw}{gives the row weights of the fourthcorner table.}
  \item{cw}{gives the column weights of the fourthcorner table.}
  \item{lw}{gives the row weights of the fourthcorner table.}
  \item{tab}{gives the fourthcorner table.}
  \item{li}{gives the coordinates of the covariates of R.}
  \item{l1}{gives the normed scores of the covariates of R.}
  \item{co}{gives the coordinates of the covariates of Q.}
  \item{c1}{gives the normed scores of the covariates of Q.}
  \item{lR}{gives the row coordinates of R.}
  \item{mR}{gives the normed row scores of R.}
  \item{lQ}{gives the row coordinates of Q.}
  \item{mQ}{gives the normed row scores of Q.}
  \item{aR}{gives projection of axis onto co-inertia axis of R.}
  \item{aR}{gives projection of axis onto co-inertia axis of Q.}
  \item{ngenes}{gives the number of analysed genes.}
}
\author{Lara Urban}
\examples{
data(Baca)
ordBaca <- ord(ExprSet = Baca, exprvar = 1, nf = 2)
ls(ordBaca)
plot(ordBaca)
}
