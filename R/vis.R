# simultaneous visualisation of stat and ord: vis

# input:
# [Stat] object of stat function;
# [Ord] object of ord function;
# [alpha] defines significance level of pvalues (default:0.05);
# [xaxis] and [yaxis] define axes of yination (default:1 and 2, resp.);
# [col] three cols: first col defines non-significant variables, second/third col
# define significant positive/negative asscoiations, respectively
# (default: "gray",transblue,transred);
# [alphatrans] defines degree of transparency of the second and third col;
# [cex] defines text size (default=1);
# [rangex, rangey] define ize of the plot (default=2);
# [...] additional features can be added

vis <- function(Stat, Ord=NULL, alpha=0.05, xaxis=1, yaxis=2,
                col=c("gray", transblue, transred),
                alphatrans=0.5, cex=1, rangex=2, rangey=2, ...) {

    if (missing(Stat)) {
        print("Only ordination will be plotted.")
    }
    else if (!inherits(Stat, "stat")) {
        stop("object of class stat is required here")
    }

    if (!inherits(Ord, "ord")) {
        stop("object of class ord is required here")
    }
    if (length(Ord$eig) < 2) stop("The fourthcorner matrix should be factorized
        into at least two eigenvectors.")

    # col: neutral, positive, negative; alphatrans: level of transparency
    transblue=rgb(0, 0, 1, alpha=alphatrans)
    transred=rgb(1, 0, 0, alpha=alphatrans)
    
    # size of text
    cex <- par("cex") * cex

    # save coordinates of rows and columns of Ord
    rowcoor <- Ord$l1[,c(xaxis,yaxis)]
    colcoor <- Ord$c1[,c(xaxis,yaxis)]

    # span the coordinate system
    plot(rangex * range(min(rowcoor[,1], colcoor[,1]), max(rowcoor[,1], colcoor[,1])),
         rangey * range(min(rowcoor[,2], colcoor[,2]), max(rowcoor[,2], colcoor[,2])),
         xlab=paste("axis", xaxis),
         ylab=paste("axis", yaxis),
         type='n',
         bty='l')

    # order rows of rowcoor and colcoor alphabetically
    rowcoor <- rowcoor[order(rownames(rowcoor)),]
    colcoor <- colcoor[order(rownames(colcoor)),]

    # determine positive and negative significant associations
    if (missing(Stat)) {
        text(rowcoor[,1], rowcoor[,2], row.names(rowcoor), cex=cex)
        text(colcoor[,1], colcoor[,2], row.names(colcoor), cex=cex)
    }
    else {
        possig <- which(Stat$adj.pvalue <= alpha & Stat$stat > 0, arr.ind=TRUE)
        negsig <- which(Stat$adj.pvalue <= alpha & Stat$stat < 0, arr.ind=TRUE)
        # save all significant associations
        sig <- list(unique(c(possig[,1], negsig[,1])),
                    unique(c(possig[,2], negsig[,2])))
        # if there are significant associations:
        # visualise segments between significant variables in cols[2,3]
        if (length(sig[[1]]) > 0) {
            if (nrow(possig) > 0) {
                segments(rowcoor[possig[,1],1], rowcoor[possig[,1],2],
                         colcoor[possig[,2],1], colcoor[possig[,2],2], lty=1, lwd=2,
                         col=col[2])
            }
            if (nrow(negsig) > 0) {
                segments(rowcoor[negsig[,1],1], rowcoor[negsig[,1],2],
                         colcoor[negsig[,2],1], colcoor[negsig[,2],2], lty=1, lwd=2,
                         col=col[3])
            }

            # visualize variables without significant associations in col[1]
            if (length(row.names(rowcoor)[-sig[[1]]]) > 0) {
                text(rowcoor[-sig[[1]],1], rowcoor[-sig[[1]],2],
                     row.names(rowcoor)[-sig[[1]]], cex=cex, col=col[1])
            }
            if (length(row.names(colcoor)[-sig[[2]]]) > 0) {
                text(colcoor[-sig[[2]],1], colcoor[-sig[[2]],2],
                     row.names(colcoor)[-sig[[2]]], cex=cex, col=col[1])
            }

            # visualize significant variables in black
            if (length(row.names(rowcoor)[sig[[1]]]) > 0) {
                text(rowcoor[sig[[1]],1], rowcoor[sig[[1]],2],
                     row.names(rowcoor)[sig[[1]]], cex=cex)
            }
            if (length(row.names(colcoor)[sig[[2]]]) > 0) {
                text(colcoor[sig[[2]],1], colcoor[sig[[2]],2],
                     row.names(colcoor)[sig[[2]]], cex=cex)
            }

            # if there are no significant associations, show all in col[1]
        } else {
            text(rowcoor[,1], rowcoor[,2], row.names(rowcoor),
                 cex=cex, col=col[1])
            text(colcoor[,1], colcoor[,2], row.names(colcoor),
                 cex=cex, col=col[1])
        }
    }
}
