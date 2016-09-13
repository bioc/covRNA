# visualisation of different features of x

# input:
# [x] object of x function;
# [feature] defines which features shall be visualised:
# "columns L","rows L", "columns R", "columns Q", "variance",
# "correlation circle R" or "correlation circle Q" (default: "variance");
# [xaxis] and [yaxis] define axes of xination (default:1 and 2, resp.);
# [cex] defines text size (default: 1);
# [rangex, rangey] define ize of the plot (default=2);
# [...] additional features can be added

plot.ord <- function(x, feature="variance", xaxis=1, yaxis=2,
                         cex=1, range=2, ...) {

  if (feature=="rows L") {
    plot(range * x$lR[, c(xaxis, yaxis)], type='n', xlab=paste("axis", xaxis),
        ylab=paste("axis", yaxis), bty='l')
    text(x$lR[,1], x$lR[,2], row.names(x$lR), cex=cex)
  }
    
  if (feature=="columns L") {
    plot(range * x$lQ[, c(xaxis, yaxis)], type='n', xlab=paste("axis", xaxis),
        ylab=paste("axis", yaxis), bty='l')
    text(x$lQ[,1], x$lQ[,2], row.names(x$lQ), cex=cex)
  }
    
  if (feature=="correlation circle R") {
    plot(range * rbind(x$aR,c(0,0)), type='n', bty='l')
    text(x$aR[,1], x$aR[,2], row.names(x$aR), cex=1)
    text(0,0, "0")
    segments(0,0,x$aR[1,1], x$aR[1,2])
    segments(0,0,x$aR[2,1], x$aR[2,2])
  }
    
  if (feature=="correlation circle Q") {
    plot(range * rbind(x$aQ,c(0,0)), type='n', bty='l')
    text(x$aQ[,1], x$aQ[,2], row.names(x$aQ), cex=1)
    text(0,0, "0")
    segments(0,0,x$aQ[1,1], x$aQ[1,2])
    segments(0,0,x$aQ[2,1], x$aQ[2,2])    
  }
    
  if (feature=="columns R") {
    plot(range * x$l1[, c(xaxis, yaxis)], type='n', xlab=paste("axis", xaxis),
         ylab=paste("axis", yaxis), bty='l')
    text(x$l1[,1], x$l1[,2], row.names(x$l1), cex=cex)
  }
    
  if (feature=="columns Q") {
    plot(range * x$c1[, c(xaxis, yaxis)], type='n', xlab=paste("axis", xaxis),
        ylab=paste("axis", yaxis), bty='l')
    text(x$c1[,1], x$c1[,2], row.names(x$c1), cex=cex)
  }
    
  if (feature=="variance") {
    barplot(x$variance, col=c(rep("black", x$nf),
                              rep("grey", length(x$variance) - x$nf)),
            names=1:length(x$variance),
            ylab="amount of explained variance",
            ylim=c(0, max(max(x$variance) + 0.1,1)))
  }
}
