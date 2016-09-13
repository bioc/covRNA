\name{plot.stat}
\alias{plot.stat}
\title{Plot the fourthcorner analysis for transcriptomic data}
\description{
The function plot produces a cross table of the gene and sample covariates of a stat object. Colours indicate positive/negative significance or absence of significance of the assciations (per default: white for non-significant, red for negative significant and red for positive significant associations).
}

\usage{
\method{plot}{stat}(x, col=c("lightgrey","deepskyblue","red"), sig=TRUE,
                    alpha=0.05, show=c("adj","non-adj"), cex=1, 
                    ynames, xnames, ytext=1, xtext=1, shiftx=0, shifty=0, ...)
}

\arguments{
  \item{x}{An object of class stat that shall be visualised as a cross table.}
  \item{col}{A vector of three colours. The first colour represents non-significant, the second positive significant, the third negative significant associations in the cross table.}
  \item{sig}{If TRUE (default), only covariates involved in at least one significant association are plotted.}
  \item{alpha}{The significance level.}
  \item{show}{'adj' or 'non-adj' indicate if adjusted or raw p-values shall be plotted, respectively.}
  \item{cex}{The magnitude of the text in the cross table.}
  \item{ynames, xnames}{Row and column names of the cross table. By default, the column names of R and Q are used, respectively.}
  \item{ytext, xtext}{Rotation of the row and column names of the cross table.}
  \item{shifty, shiftx}{Shift of the row and column names to the right or to the left.}
  \item{...}{More plotting parameters can be added.}
}
\value{Plot of fourthcorner analysis.}
\author{Lara Urban}
\examples{
statBaca <- stat(Baca, nrcor = 2)
plot(statBaca)
}


