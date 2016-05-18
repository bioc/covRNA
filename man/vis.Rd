\name{vis}
\alias{vis}
\title{Simultaneous visualisation of transcriptomic data by combining fourthcorner analysis and RLQ}
\description{
The vis function simultaneously visualizes the results of the functions stat and ord. Firstly, all covariates of R and Q are visualized by ordination in one plot; covariates involved in at least one significant association are shown in black, other covariates are shown in gray. Then, all covariates that are significantly associated according to stat are connected by lines which color represents the character of their significance.
}
\usage{
vis(Stat, Ord=NULL, alpha=0.05, xaxis=1, yaxis=2, col=c("gray", transblue, transred),
    alphatrans=0.5, cex=1, rangex=2, rangey=2, ...)
}
\arguments{
  \item{Stat}{An object of class stat.}
  \item{Ord}{An object of class ord. The objects stat and ord should have the same value ngenes.}
  \item{alpha}{The significance level.}
  \item{xaxis, yaxis}{Define which axes of ordination shall be shown by x- and y-axis, respectively.}
  \item{col}{A vector of three colors. The first color represents non-significant variables, the second positive significant, the third negative significant associations.}
  \item{alphatrans}{Defines degree of transparency of the second and third color.}
  \item{cex}{The magnitude of the text in the ordination.}
  \item{rangex, rangey}{The range of the x axis and y axis can be extended or reduced, e.g. for the case that not all covariates are visible in the default setting.}
  \item{...}{More plotting parameters can be added.}
}
\value{Plot of fourthcorner analysis and RLQ. }
\author{Lara Urban}
\examples{
data(Baca)
statBaca <- stat(Baca, nrcor = 2)
ordBaca <- ord(Baca)
vis(Stat = statBaca, Ord = ordBaca)
vis(Ord = ordBaca)
}
