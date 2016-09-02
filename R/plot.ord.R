# visualisation of different features of x

# input:
# [x] object of x function;
# [feature] defines which features shall be visualised:
# "columns L","rows L", "columns R", "columns Q", "variance",
# "correlation circle R" or "correlation circle Q" (default: "variance");
# [xaxis] and [yaxis] define axes of xination (default:1 and 2, resp.);
# [labelsize] defines size of the label (default: 1.25);
# [boxes] if boxes shall be drawn around labels (default: TRUE);
# [...] additional features can be added

plot.ord <- function(x, feature="variance", xaxis=1, yaxis=2,
                         labelsize=1.25, boxes=TRUE, ...) {

  if (feature=="rows L") {
    plot(x$lR[, c(xaxis, yaxis)], type='n', xlab=paste("axis", xaxis),
        ylab=paste("axis", yaxis), bty='l')
    text(x$lR[,1], x$lR[,2],
        row.names(x$lR), cex=labelsize)
    #s.label(x$lR[, c(xaxis, yaxis)], clabel=labelsize, boxes=boxes)
  }
  if (feature=="columns L") {
    s.label(x$lQ[, c(xaxis, yaxis)], clabel=labelsize, boxes=boxes)
  }
  if (feature=="correlation circle R") {
    s.corcircle(x$aR, xaxis, yaxis, clabel=labelsize)
  }
  if (feature=="correlation circle Q") {
    s.corcircle(x$aQ, xaxis, yaxis, clabel=labelsize)
  }
  if (feature=="columns R") {
    s.arrow(x$l1, xax=xaxis, yax=yaxis, clabel=labelsize, boxes=boxes)
  }
  if (feature=="columns Q") {
    s.arrow(x$c1, xax=xaxis, yax=yaxis, clabel=labelsize, boxes=boxes)
  }
  if (feature=="variance") {
    barplot(x$variance, col=c(rep("black", x$nf),
                                    rep("grey", length(x$variance) - x$nf)),
            names=1:length(x$variance),
            ylab="amount of explained variance",
            ylim=c(0, max(max(x$variance) + 0.1,1)))
  }
}
