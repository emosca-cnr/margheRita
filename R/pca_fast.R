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

pca_fast <- function(mRList, dirout, col_by="class", method="svd", scaling=c("none", "pareto", "vector", "uv"), center=TRUE, include_QC=TRUE, top=Inf, nPcs, write_output=FALSE, ...) {


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

  #Graphic pairs

  pairplot= paste(dirout, "/Pairs.png", sep = "")
  grDevices::png(
    pairplot,
    width = 6,
    height = 6,
    units = "in",
    res= 300
  )

  col_factor <- as.factor(X_ann[, col_by])
  col_pal <- rainbow(length(levels(col_factor)))
  pairs(mRList$pca@scores[,1:nPcs],labels=paste(colnames(mRList$pca@scores),"(",round(mRList$pca@R2*100,3), "%)"),
        main = "Pairs",
        col = col_pal[as.numeric(col_factor)],
        pch = 19
  )
  par(xpd = TRUE)
  legend("topleft", legend = levels(col_factor), col = col_pal, pch=8, cex=0.5,ncol=6)

  grDevices::dev.off()

  if (write_output){
    scoretable= paste(dirout, "/Scoretable.csv", sep ="")
    utils::write.csv(mRList$pca@scores, scoretable)
    loadingtable = paste(dirout ,"/Loadingtable.csv", sep= "")
    utils::write.csv(mRList$pca@loadings, loadingtable)
    variance= paste(dirout, "/Variance.csv", sep = "")
    utils::write.csv(mRList$pca@R2cum,variance)

  }

  return(mRList)

}
