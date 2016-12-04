# statistical analysis of covRNA: stat

# input:
# [ExprSet] ExpressionSet: assayData will be called L (non-negative
# values), phenoData Q and featureData R (default);
# or: [R, L, Q] single dataframes;
# [npermut] number of permutations;
# [padjust] multiple hypothesis testing adjustment method;
# [nrcor] number of cores that shall be used (default:all available-1);
# [exprvar] for gene filtering: may vary between 0 and 1 (=default)
# & indicates fraction of most variably expressed genes to be analysed;

# output: list of class stat:
# [pvalue], [adj.pvalue], [stat] matrices of adjusted, non-adjusted
# p-values and statistical values;
# [ngenes] number of analyzed genes;
# [adjust.method] multiple hypothesis testing adjustment method
# (default: Benjamini and Hochberg)
# [npermut] number of permutations (default: 9999)
# [call] shows called function

# statistics:
# quant-quant: correlation coefficient of normalised variables
# qual-qual: Chi-square contingency table of binary variables
# qual-quant: correlation coefficient

stat <- function(ExprSet, R=NULL, L=NULL, Q=NULL, npermut=9999, padjust="BH", nrcor=detectCores()-1, exprvar=1) {

  #####

  # subfunctions:

  # 1. "disj" returns the disjunctive table of a factor vector
  disj <- function(fac) {
    # initialise matrix with levels of the factor as columns
    mat <- matrix(0, nrow = length(fac), ncol = nlevels(fac))
    # set respective values to 1, the rest remains 0
    mat[(1:length(fac)) + length(fac) * (unclass(fac)-1)] <- 1
    # give row- and columnames to the matrix
    dimnames(mat) <- list(names(fac), as.character(levels(fac)))
    # return the matrix as dataframe
    return(data.frame(mat, check.names = FALSE))
  }

  # 2. "normalise" normalises R/Q by the row/column weights of L
  normalise <- function (mtx, weight, sumL) {
    # mtx as matrix R or Q to normalise, weight as row/column weight of
    # matrix L and sumL as sum of entire matrix L
    # row-wise mean of mtx
    mn <- 0
    mn <- sum(mtx * (weight/sumL))
    # row-wise variance of mtx
    vr <- 0
    vr <- sum((weight/sumL) * (mtx-mn)^2)
    # set variances <= 0 to 1 to avoid effect
    if (vr<=0) vr <- 1
    # calculate standard deviation
    vr <- sqrt(vr)
    # normalise each column of mtx
    mtx <- (mtx-mn)/vr
    return(mtx)
  }

  # 3. "calstat" permutes matrix L row- or columnwise once and counts how
  # often empirical statistical values >= observed statistical value
  calstat <- function(R, Rq, L, Q, Qq, o, permut, ind) {
    # permut="r" indicates row-wise permutation
    # permut="c" indicates column-wise permutation
    if (permut=="r") {Lpermut <- L[sample(1:nrow(L)),]}
    if (permut=="c") {Lpermut <- L[,sample(1:ncol(L))]}
    # calculate permuted statistic of quantitative and categorical
    newobsc <- t(R) %*% (Lpermut %*% Q)
    newobsq <- t(Rq) %*% (Lpermut %*% Qq)
    # dependent on index, decide which observed value shall be chosen
    newobs <- matrix(sapply(1:prod(dim(ind)), function(x)
                        if (as.numeric(ind)[x]==4) as.numeric(newobsc)[x]
                        else as.numeric(newobsq)[x]),
                     nrow=nrow(ind),
                     ncol=ncol(ind))
    # count how often permuted statistic is larger than observed one
    return(newobs>=o)
  }


  #####

  # examine and transform data

  # if ExpressionSet is missing: use single dataframes as input
    if (missing(ExprSet)) {
        if (missing(R)) {
            R = as.data.frame(diag(x=1, nrow=nrow(L), ncol=nrow(L)))
            warning("R is missing and is replaced by an identity matrix.")
        }
        if (missing(Q)) {
            Q = as.data.frame(diag(x=1, nrow=ncol(L), ncol=ncol(L)))
            warning("R is missing and is replaced by an identity matrix.")
        }
        if (missing(L)) {
            stop("L is missing.")
        }

        if (is.matrix(L)) matL <- L else
            if (is.data.frame(L)) matL <- as.matrix(L) else
                stop("L has to be a dataframe or a matrix")

            if (is.matrix(R)) {
                tabR <- as.data.frame(R)
                warning("R is transformed to a dataframe by as.data.frame().")
            } else if (is.data.frame(R)) tabR <- R else
                stop("R has to be a dataframe or a matrix")

            if (is.matrix(Q)) {
                tabQ <- as.data.frame(Q)
                warning("Q is transformed to a dataframe by as.data.frame().")
            } else if (is.data.frame(Q)) tabQ <- Q else
                stop("Q has to be a dataframe or a matrix")
    } else {
        # else: transform ExpressionSet to matrix/dataframes
    matL <- exprs(ExprSet)
    tabR <- fData(ExprSet)
    tabQ <- pData(ExprSet)
  }

  # control values of the dataframes
  if (any(is.na(tabR)))
    stop("NA in R")
  if (any(is.na(tabQ)))
    stop("NA in Q")
  if (any(is.na(matL)))
    stop("NA in L")
  if (any(matL<0))
    stop("negative values in L")

  # control features of the dataframes
  if (nrow(tabR) != nrow(matL))
    stop("row number of R does not match row number of L")
  if (nrow(tabQ) != ncol(matL))
    stop("row number of Q does not match column number of L")

  # filter gene by gene expression variance according to exprvar
  if (exprvar != 1) {
    # determine number of genes by given fraction
    ngenes <- trunc(exprvar * nrow(matL))
    # filter L to retain the most variably expressed genes
    matL <- matL[order(rowVars(matL), decreasing=TRUE)[1:ngenes],]
    # just use gene covariates in R which contain at least one annotation after
    # the unsupervised filtering
    whichR <- unique(which(colSums(tabR)==0))
    tabR <- tabR[,-whichR]
  }

  # transform R and Q to matrices and order columns alphabetically
  # R:
  # init matrix
  matR <- matrix(0, nrow(tabR), 1)
  # init colnames of matrix
  namesR <- "init"
  # init indexR which will assign 1 to quantitative and 2 to
  # categorical variables
  indexR <- "init"
  for (i in 1:ncol(tabR)) {
    # if variable is quantitative: add it to new matrix
    if (is.numeric(tabR[,i])) {
      matR <- cbind(matR, tabR[,i])
      namesR <- c(namesR, names(tabR)[i])
      indexR <- c(indexR, 1)
    }
    # if variable is categorical: add disjunctive table to matrix
    else if (is.factor(tabR[,i])) {
      disjtabR <- disj(tabR[,i])
      matR <- cbind(matR, disjtabR)
      # colnames of the variables will contain variable and level name
      namesR <- c(namesR, paste(names(tabR)[i], ".", names(disjtabR)))
      indexR <- c(indexR, rep(2, ncol(disjtabR)))
    }
  }
  # delete the initialisations
  matR <- as.matrix(matR[,-1])
  colnames(matR) <- namesR[-1]
  rm(namesR)
  indexR <- indexR[-1]
  # order the variables alphabetically
  if (dim(matR)[2] !=1) {
      indexR <- indexR[order(colnames(matR))]
      matR <- matR[, order(colnames(matR))]
      }
  
  


  # Q: do the same with matrix Q (see matrix R)
  matQ <- matrix(0, nrow(tabQ), 1)
  namesQ <- "init"
  indexQ <- "init"
  for (i in 1:ncol(tabQ)) {
    if (is.numeric(tabQ[,i])) {
      matQ <- cbind(matQ, tabQ[,i])
      namesQ <- c(namesQ, names(tabQ)[i])
      indexQ <- c(indexQ, 1)
    }
  else if (is.factor(tabQ[,i])) {
    disjtabQ <- disj(tabQ[,i])
    matQ <- cbind(matQ, disjtabQ)
    namesQ <- c(namesQ, paste(names(tabQ)[i], ".", names(disjtabQ)))
    indexQ <- c(indexQ, rep(2, ncol(disjtabQ)))
    }
  }
  matQ <- as.matrix(matQ[,-1])
  colnames(matQ) <- namesQ[-1]
  rm(namesQ)
  indexQ <- indexQ[-1]
  if (dim(matQ)[2] !=1) {
      indexQ <- indexQ[order(colnames(matQ))]
      matQ <- matQ[, order(colnames(matQ))]
      }
  

  # calculate index as sum of indexR and indexQ; only index==4 is interesting,
  # since the both variables are categorical (important for
  # calculation of statistical values)
  index <- matrix(as.numeric(rep(indexR, length(indexQ))),
                  nrow=length(indexR), ncol=length(indexQ))
  index <- matrix(apply(index, 1, function(x) (x + as.numeric(indexQ))),
                    nrow = length(indexR), ncol = length(indexQ), byrow=TRUE)

  # normalise quantitative variables of R and Q per column and save these
  # in matRq and matQq
  matRq <- apply(matR, 2, function(x) normalise(x, rowSums(matL), sum(matL)))
  matQq <- apply(matQ, 2, function(x) normalise(x, colSums(matL), sum(matL)))

  #####

  # results in the list result

  result <- list()

  # calculate observed statistic for quantitative and categorical matrices
  obsc <- t(matR)%*%(matL%*%matQ)
  obsq <- t(matRq)%*%(matL%*%matQq)

  # if both variables are categorical (index==4), cell of obsc is taken;
  # if at least one of the variables is quantitative: take obsq
  obs <- matrix(sapply(1:prod(dim(index)), function(x)
      if (as.numeric(index)[x] == 4) as.numeric(obsc)[x] else
          as.numeric(obsq)[x]), nrow=nrow(index), ncol=ncol(index))

  # create as many clusters as number of available cores-1
  clus <- makeCluster(nrcor)

  # export needed data to clusters
  clusterExport(clus,list("matL","matR","matQ","matRq","matQq","index",
                          "obs","calstat"),
                envir=environment())

  # use parLapply to count how often empirical statistical values >= obs
  # in the case of row-wise permutation
  countr <- parLapply(clus, rep("r",npermut), function(x)
                        calstat(R=matR, Rq=matRq, L=matL, Q=matQ, Qq=matQq,
                        o=obs, permut=x, ind=index))
  countr <- Reduce(function(x,y){x+y}, countr)
  # in the case of column-wise permutation
  countc <- parLapply(clus, rep("c",npermut), function(x)
                        calstat(R=matR, Rq=matRq, L=matL, Q=matQ, Qq=matQq,
                        o=obs, permut=x, ind=index))
  countc <- Reduce(function(x,y){x+y}, countc)

  # close the clusters
  stopCluster(clus)

  # save the statistical value dependent on the index
  # if both variables are categorical (index==4), the N statistic amounts to
  # the already calculated statistic
  # if at least one variables is quantitative, the correlation coefficient is
  # computed by dividing the statistic by the sum over matrix L
  result$stat <- matrix(sapply(1:prod(dim(index)), function(x)
                            if (as.numeric(index)[x] == 4) {as.numeric(obs)[x]}
                            else {as.numeric(obs)[x]/sum(matL)}),
                        nrow=nrow(index), ncol=ncol(index),
                        dimnames=list(colnames(matR), colnames(matQ)))

  # calculate p-value of row-wise and column-wise permutation
  # independently from each other (equation see Vignette)
  pvaluer <- 2 * pmin((countr+1)/(npermut+1), (1-(countr+1)/(npermut+1)))
  pvaluec <- 2 * pmin((countc+1)/(npermut+1), (1-(countc+1)/(npermut+1)))

  # special case if variance of at least one variable equals zero:
  # all pvalues of this variable are 1
  pvaluer[which(rowVars(t(matR))==0),] <- 1
  pvaluer[,which(rowVars(t(matQ))==0)] <- 1
  pvaluec[which(rowVars(t(matR))==0),] <- 1
  pvaluec[,which(rowVars(t(matQ))==0)] <- 1

  # adjust p-values and choose maximal pvalue
  adj.pvaluer <- as.numeric(p.adjust(pvaluer, method=padjust))
  adj.pvaluec <- as.numeric(p.adjust(pvaluec, method=padjust))
  result$adj.pvalue <- matrix(pmax(adj.pvaluer, adj.pvaluec),
                              ncol=ncol(result$stat),
                              nrow=nrow(result$stat),
                              dimnames=list(rownames(result$stat),
                                            colnames(result$stat)))

  # investigate if the maximum p-value is taken from row- or column-wise
  # permutation and choose respective non-adjusted p-value
  is.r=(result$adj.pvalue == adj.pvaluer)
  pvalue <- rep(0, length(is.r))
  for (i in 1:length(is.r)) {
    if (is.r[i] == TRUE) {pvalue[i] <- pvaluer[i]}
    else {pvalue[i] <- pvaluec[i]}
  }
  
  # statistical tests
  stattest <- matrix('corr coeff', dimnames = list(colnames(matR), colnames(matQ)), 
                     nrow = length(indexR), ncol=length(indexQ))
  stattest[index==4] <- 'Chi-square related test'
  
  # assign p-values
  result$pvalue <- matrix(pvalue, ncol=ncol(result$stat),
                          nrow=nrow(result$stat),
                          dimnames=list(rownames(result$stat),
                                        colnames(result$stat)))

  # add call of the function to result
  result$call <- match.call()
  result$stattest <- stattest
  # add number of analysed genes to result
  if (exprvar != 1) {result$ngenes <- ngenes} else result$ngenes <- "all"
  # add number of permutations to result
  result$npermut <- npermut
  # add multiple hypothesis adjustment method to result
  result$adjust.method <- padjust
  class(result) <- "stat"
  return(result)
}




