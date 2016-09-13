# visualisation of the stat results as table of covariate combinations

# input:
# [x] object of stat function;
# [col] three colours: first colour defines non-significant associations,
# second/third colour define significant positive/negative asscoiations
# (default: "lightgrey","deepskyblue","red");
# [sig] if TRUE, only covariates involved in at least one significant
# association are shown (default: TRUE)
# [alpha] defines significance level of pvalues (default:0.05);
# [show] defines which p-values shall be chosen (default:adj);
# [cex] defines text size of table;
# [ynames] and [xnames] define row and column names of the table (optional);
# [ytext] and [xtext] define how text shall be rotated (default: 1 for
# horizontal; 3 for vertical);
# [shiftx] and [shifty] define how much the text should be shifted to the
# right/top (if>0) or to the left/bottom (if<0);
# [...] additional features can be added

plot.stat <- function(x, col=c("lightgrey","deepskyblue","red"), sig=TRUE,
                      alpha=0.05, show=c("adj","non-adj"), cex=1,
                      ynames, xnames, ytext=1, xtext=1, shiftx=0, shifty=0,
                      ...) {

  # save statistical values as matrix
  stat <- x$stat
  statnum <- as.numeric(stat)

  # save either adjusted or non-adjusted p-values as matrix (default:adj)
  show <- match.arg(show)
  if (show=="non-adj") {
    statpvalm <- x$pvalue
  } else {
    statpvalm <- x$adj.pvalue 
  }
  
  if (sig == TRUE) {
    rows <- unique(which(statpvalm<=alpha, arr.ind = TRUE)[,1])
    columns <- unique(which(statpvalm<=alpha, arr.ind = TRUE)[,2])
    statpvalm <- statpvalm[rows, columns]
    statpval <- as.numeric(statpvalm)
    stat <- stat[rows, columns]
    statnum <- as.numeric(stat)
  }

  # choose significant associations
  statsel <- which(statpval <= alpha)

  # init vector des with 1; if des==1: no significant association
  des <- rep(1, length(statnum))

  # set des=2 for postive significant associations and des=3 for
  # negative significant associations
  des[statsel] <- ifelse(statnum[statsel] > 0,2,3)

  # define parameters for plotting a table
  xx <- 1:ncol(stat)
  yy <- nrow(stat):1
  wx <- range(xx)
  wy <- range(yy)
  dx <- diff(wx)/(length(xx)-1)/2
  dy <- diff(wy)/(length(yy)-1)/2
  xlim <- wx + c(-dx, dx)
  ylim <- wy +c (-dy, dy)

  # set graphical parameters
  strx <- max(strwidth(paste(" ", rownames(stat), " ", sep =""),
                         units="inches", cex=cex))
  stry <- max(strwidth(paste(" ", colnames(stat), " ", sep =""),
                         units = "inches", cex=cex))
  par(mai=c(0.2, strx, stry, 0.1))

  # plot an empty table
  plot.default(0, 0, type="n", xlab="", ylab="",
               xaxt="n", yaxt="n", xlim=xlim, ylim=ylim,
               xaxs="i", yaxs="i", frame.plot=FALSE)

  # add text to the rows and columns of the table
  # if no names are given as input, take row- and columnnames of stat
  if (missing(ynames)) ynames <- rownames(stat)
  if (missing(xnames)) xnames <- colnames(stat)
  # plot names next to table
  ynew <- seq(wy[1], wy[2], le=length(yy))[rank(yy)]
  mtext(at=ynew+shifty, side=2, text=paste(ynames," ",sep=""), adj=1,
        cex=cex, las=ytext)
  xnew <- seq(wx[1], wx[2], le=length(xx))[rank(xx)]
  mtext(at=xnew+shiftx, side=3, text=paste(xnames," ",sep=""), adj=0,
        cex=cex, las=xtext)

  # add coloured rectangles to the table
  xtot <- xx[col(stat)]
  ytot <- yy[row(stat)]
  xdelta <- 0.9*(max(xx)-min(xx))/(length(xx)-1)/2
  ydelta <- 0.9*(max(yy)-min(yy))/(length(yy)-1)/2
  rect(xtot-xdelta, ytot-ydelta, xtot+xdelta, ytot+ydelta, col=col[1:3][des],
       border="black")
}



