#' PCA
#' Function that allows to perform PCA. Graphs of Scree plot and pairs of the first 10 components were plotted. Score and loading plot, choosing pcx and pcy, were plotted launching function plot2DPCA after pca function.
#' @param mRList mRList object
#' @param col_by (default class)
#' @param method see pcaMethods::pca()
#' @param center whether to center the data or not
#' @param scaling see pcaMethods::prep()
#' @param include_QC (default TRUE)
#' @return  mRList object with "pca" element
#' @export
#' @importFrom graphics plot barplot pairs
#' @importFrom grDevices dev.off png rainbow
#' @importFrom utils write.csv
#' @importFrom stats var prcomp
#' @importFrom pcaMethods pca
#' @param top only the top most varying features will be used
#' @param dirout output directory
#' @param nPcs number of principal components
#' @param ... further argments to pcaMethods::pca

pca_fast <- function(mRList, dirout, col_by="class", method="svd", scaling=c("none", "pareto", "vector", "uv"), center=TRUE, include_QC=TRUE, top=Inf, nPcs=5, ...) {
  
  
  scaling <- match.arg(scaling)
  
  dir.create(dirout)
  
  X <- mRList$data
  X_ann <- mRList$sample_ann
  
  #include QC
  if(include_QC){
    cat("Including QC\n")
    X <- cbind(X, mRList$QC)
    X_ann<- rbind(X_ann, mRList$QC_ann)
  }
  
  if(nrow(X) > top){
    cat("using only top", top, "metabolites by variance\n")
    idx <- order(-apply(X, 1, var))[1:top]
    X <- X[idx, ]
  }
  
 
  #it is the same as PCA for euclidean distance
  
  mRList$pca <- pcaMethods::pca(t(X), method = method, nPcs = nPcs, scale = scaling, center = center, ...)
 
  
  #Graphic screeplot
  scree = paste(dirout, "/Screeplot.png", sep = "")
  grDevices::png(
    scree,
    width = 20,
    height = 20,
    units = "cm",
    res = 300
  )
  graphics::barplot(
    mRList$pca@R2,
    xlab = "Principal Components",
    ylab = "Proportion of Variance explained",
    main = "Screeplot"
  )
  grDevices::dev.off()
  
  
  return(mRList)
  
}
