#' Principal Component Analysis
#' @param mRList mRList object
#' @param col_by level of column class to use for assigning colors
#' @param method see pcaMethods::pca()
#' @param center whether to center the data or not
#' @param scaling see pcaMethods::prep()
#' @param include_QC (default TRUE)
#' @param nPcs number of principal components (PCs)
#' @importFrom graphics par plot barplot pairs
#' @importFrom grDevices dev.off png rainbow
#' @importFrom pcaMethods pca
#' @param top only the top most varying features will be used
#' @param dirout output directory, if not null in this directory will be saved csv files with scores, loadings and R2.
#' @param nPcs number of principal components
#' @param ... further argments to `pcaMethods::pca()`
#' @return  mRList object with "pca" element. 
#' @export

mR_pca <- function(mRList=NULL, dirout=NULL, col_by="class", method="svd", scaling=c("none", "pareto", "vector", "uv"), center=TRUE, include_QC=TRUE, top=Inf, nPcs=2, ...) {


  scaling <- match.arg(scaling)

  if (!is.null(dirout)) {
    dir.create(dirout, showWarnings = F)
  }else{
    dirout <- getwd()
  }
  
  X <- mRList$data
  X_ann <- mRList$sample_ann

  #include QC
  if(include_QC){
    cat("Including QC\n")
    X <- cbind(X, mRList$QC)
    X_ann <- rbind(X_ann, mRList$QC_ann)
  }

  if(nrow(X) > top){
    cat("using only top", top, "metabolites by variance\n")
    idx <- order(-apply(X, 1, var))[1:top]
    X <- X[idx, ]
  }


  #it is the same as PCA for euclidean distance

  mRList$pca <- pca(t(X), method = method, nPcs = nPcs, scale = scaling, center = center, ...)

  write.csv(mRList$pca@scores, file.path(dirout, "scores.csv"))
  write.csv(mRList$pca@loadings, file.path(dirout, "loadings.csv"))
  write.csv(mRList$pca@R2cum, file.path(dirout, "variance.csv"))

  return(mRList)

}
