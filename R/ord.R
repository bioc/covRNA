# rlq analysis of covaRNA: ord

# input:
# [ExprSet] ExpressionSet: assayData will be called L (non-negative
# values), phenoData Q and featureData R (default);
# or: [R, L, Q] single dataframes;
# [exprvar] for gene filtering: may vary between 0 and 1 (=default)
# & indicates fraction of most variably expressed genes to be analysed;
# [nf] number of axes to be considered by ordination

# output: list of class ord:
# all features of an rlq object (see ade4)
# [ngenes] number of analyzed genes;
# [call] shows called function
# [variance] variance explained by the considered axes in [%]

ord <- function(ExprSet, R=NULL, L=NULL, Q=NULL, exprvar=1, nf=2) {

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
    stop("row number of R does not suit row number of L")
  if (nrow(tabQ) != ncol(matL))
    stop("row number of Q does not suit column number of L")

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

  # get variable type of R and Q
  varR <- vector(mode="logical", length=dim(tabR)[2])
  for (i in 1:dim(tabR)[2]) {
      varR[i] <- is.numeric(tabR[,i])
  }
  
  varQ <- vector(mode="logical", length=dim(tabQ)[2])
  for (i in 1:dim(tabQ)[2]) {
      varQ[i] <- is.numeric(tabQ[,i])
  }

  # visualize the three tables with respect to variable type of R/Q

  # matL is visualised by CA
  VL <- dudi.coa(matL, scannf=FALSE)

  # tabR and tabQ are visulaised by PCA, Hillsmith or MCA
  if (all(varR)) VR <- dudi.pca(tabR, row.w=VL$lw, scannf=FALSE) else 
  if (any(varR)) VR <- dudi.hillsmith(tabR, row.w=VL$lw, scannf=FALSE) else 
  VR <- dudi.acm(tabR, row.w=VL$lw, scannf=FALSE)

  if (all(varQ)) VQ <- dudi.pca(tabQ, row.w=VL$cw, scannf=FALSE) else 
  if (any(varQ)) VQ <- dudi.hillsmith(tabQ, row.w=VL$cw, scannf=FALSE) else 
  VQ <-dudi.acm(tabQ, row.w=VL$cw, scannf=FALSE)

  # RLQ is based on the three singular matrix ordinations
  result <- rlq(VR, VL, VQ, scannf=FALSE, nf=2)

  # amount of variance explained by the axes
  result$variance <- (100*result$eig/sum(result$eig))
  # add number of analysed genes to result
  if (exprvar != 1) {result$ngenes <- ngenes} else result$ngenes <- "all"
  # add call of the function to result
  result$call <- match.call()
  class(result) <- "ord"
  return(result)
}

