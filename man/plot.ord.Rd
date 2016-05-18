\name{plot.ord}
\alias{plot.ord}
\title{Plot RLQ for transcriptomic data}
\description{
The function plot can visualise different features of an ord object by adjusting the argument "feature". By default, a barplot of the variance explained by the axes of the RLQ is plotted (see arguments).
}
\usage{
\method{plot}{ord}(x, feature="variance", xaxis=1, yaxis=2, labelsize=1.25, boxes=TRUE, ...)
}
\arguments{
  \item{x}{An object of class ord that shall be visualised by ordination.}
  \item{feature}{Defines which features of the object shall be visualised: "columns L","rows L", "columns R" and "columns Q" visualise the respective variables as oridnation, "variance" shows a barplot of the variance explained by the axes, "correlation circle R" and "correlation circle Q" visualise the projection of the original space into the ordination space.}
  \item{xaxis, yaxis}{Define which axes of ordination shall be shown by x- and y-axis, respectively.}
  \item{labelsize}{Defines size of label of covariates.}
  \item{boxes}{Defines if boxes shall be drawn around labels of covariates.}
  \item{...}{More plotting parameters can be added.}
}
\value{Plot of RLQ.}
\author{Lara Urban}
\examples{
ordBaca <- ord(Baca)
plot(ordBaca)
}

