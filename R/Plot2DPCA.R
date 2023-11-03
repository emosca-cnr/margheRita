#' Plot score and loading plot
#' @description Plot 2D PCA score and loading plots choosing the principal components.
#' @param mRList mRList object
#' @param pcx Specify the principal component on the x-axis
#' @param pcy Specify the principal component on the y-axis
#' @param dirout output directory
#' @param col_by (default class)
#' @param include_QC (default TRUE)
#' @export
#' @importFrom pcaMethods scores
#' @importFrom grDevices png rainbow dev.off

Plot2DPCA <- function(mRList=NULL, pcx=1, pcy=2, dirout = "./", col_by="class", include_QC=FALSE){
  
  dir.create(dirout)
  
  xlabel = paste("PC", pcx, "(",round(mRList$pca@R2[pcx]*100,3), "%)")
  ylabel = paste("PC", pcy, "(",round(mRList$pca@R2[pcy]*100,3), "%)")
  xy <- scores(mRList$pca)
  xy <- xy[, c(pcx, pcy)]
  
  #include QC
  X_ann <- mRList$sample_ann
  if(include_QC){
    cat("Including QC\n")
    X_ann<- rbind(X_ann, mRList$QC_ann)
  }
  
  #ensure the same order between annotation and data
  X_ann <- X_ann[match(rownames(xy), X_ann[, 1]), ]
  
  col_factor <- as.factor(X_ann[, col_by])
  col_pal <- rainbow(length(levels(col_factor)))
  
  png(file.path(dirout, "scores.png"), width = 200, height = 200, units = "mm", res = 300)
  
  plot(
    xy[, 1], xy[, 2],
    xlab=xlabel, ylab=ylabel,
    main = "Score plot",
    col = col_pal[as.numeric(col_factor)],
    pch = 19
  )
  legend("bottomright", legend = levels(col_factor), col = col_pal, pch=16, cex=1)
  
  dev.off()
  
  xlabeL = paste("Loading", pcx)
  ylabeL = paste("Loading", pcy)
  xyL <- as.data.frame(mRList$pca@loadings)
  xyL <- xyL[, c(pcx, pcy)]
  

  png(file.path(dirout, "loadings.png"), width = 200, height = 200, units = "mm", res = 300)
  
  plot(
    xyL[, 1], xyL[, 2],
    xlab=xlabeL, ylab=ylabeL,
    main = "Loading plot",
    pch = 19
  )
  #text(xyL[,1], xyL[,2], labels=rownames(xyL), cex=0.8, pos=3)
  
  dev.off()
}
