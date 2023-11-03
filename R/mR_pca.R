#' PCA
#' @description This function allows to perform PCA.
#' @param mRList mRList object
#' @param col_by (default class)
#' @param method see pcaMethods::pca()
#' @param center whether to center the data or not
#' @param scaling see pcaMethods::prep()
#' @param include_QC (default TRUE)
#' @param nPcs number of principal components (PCs)
#' @param write_output if score and loading table have to be saved (default=FALSE)
#' @return  mRList object with "pca" element. Images of screeplot and pairs, table of score,loading and variance can be saved if write_output=TRUE.
#' @export
#' @importFrom graphics plot barplot
#' @importFrom grDevices dev.off png rainbow
#' @importFrom pcaMethods pca
#' @param top only the top most varying features will be used
#' @param dirout output directory
#' @param nPcs number of principal components
#' @param ... further argments to pcaMethods::pca
#' @examples
#' ##library(dataset.margheRita)
#' ##dataset(norm_pos)
#' pca_fast(mRList, col_by="class", scaling="pareto",nPcs=5,write_output=FALSE)

mR_pca <- function(mRList=NULL, dirout="./", col_by="class", method="svd", scaling=c("none", "pareto", "vector", "uv"), center=TRUE, include_QC=TRUE, top=Inf, nPcs=2, ...) {
  
  
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
  png(file.path(dirout, "scree.png"), width = 20, height = 20, units = "cm", res = 300)
  
  barplot(
    mRList$pca@R2,
    xlab = "Principal Components",
    ylab = "Proportion of Variance explained",
    main = "Screeplot"
  )
  dev.off()
  
  #Graphic pairs
  
  grDevices::png(
    file.path(dirout, "pairs.png"),
    width = 20,
    height = 20,
    units = "cm",
    res= 300
  )
  
  col_factor <- as.factor(X_ann[, col_by])
  col_pal <- rainbow(length(levels(col_factor)))
  pairs(mRList$pca@scores[,1:nPcs],labels=paste(colnames(mRList$pca@scores),"(",round(mRList$pca@R2*100,3), "%)"),
        main = "Pairs",
        col = col_pal[as.numeric(col_factor)],
        pch = 19,
        oma= c(3,3,3,15)
  )
  
  par(xpd = TRUE)
  legend("bottomright", legend = levels(col_factor), col = col_pal, pch=8, cex=0.8,ncol=1)
  
  dev.off()
  
  write.csv(mRList$pca@scores, file.path(dirout, "scores.txt"))
  write.csv(mRList$pca@loadings, file.path(dirout, "loadings.txt"))
  write.csv(mRList$pca@R2cum, file.path(dirout, "variance.txt"))
  
  return(mRList)
  
}
