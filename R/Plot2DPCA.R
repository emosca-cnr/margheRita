#' Create 2D PCA score and loading plots
#' Create 2D PCA score and loading plots choosing the principal components.
#' @param mRList mRList object
#' @param pcx Specify the principal component on the x-axis
#' @param pcy Specify the principal component on the y-axis
#' @param dirout output directory
#' @param col_by (default class)
#' @param include_QC (default TRUE)
#' @export
#' @importFrom graphics plot legend
#' @importFrom grDevices png rainbow

Plot2DPCA <- function(mRList, pcx, pcy, dirout = "./", col_by="class", include_QC=TRUE){
  
  dir.create(dirout)
  
  xlabel = paste("PC",pcx, "(", mRList$pca$variance[pcx],1, "%)")
  ylabel = paste("PC",pcy, "(", mRList$pca$variance[pcy],1, "%)")
  pc1 = mRList$pca$x[, pcx]
  pc2 = mRList$pca$x[, pcy]
  
  X <- mRList$data
  X_ann <- mRList$sample_ann
  
  #include QC
  if(include_QC){
    cat("Including QC\n")
    X <- cbind(X, mRList$QC)
    X_ann<- rbind(X_ann, mRList$QC_ann)
  }
  
  
  scoreplot = paste(dirout, "/Scoreplot.png", sep = "")
  grDevices::png(
    scoreplot,
    width = 8,
    height = 8,
    units = "in",
    res = 300
  )
  
  col_factor <- as.factor(X_ann[, col_by])
  col_pal <- rainbow(length(levels(col_factor)))
  graphics::plot(
    pc1,pc2,
    xlab=xlabel,ylab=ylabel,
    main = "Score plot",
    col = col_pal[as.numeric(col_factor)],
    pch = 19
  )
  legend("bottomright", legend = levels(col_factor), col = col_pal, pch=16, cex=0.5)
  
  grDevices::dev.off()
  
  
  xlabeL = paste("Loading",pcx)
  ylabeL = paste("Loading",pcy)
  pc1L = mRList$pca$rotation[, pcx]
  pc2L = mRList$pca$rotation[, pcy]
  
  loadingplot = paste(dirout, "/Loadingplot.png", sep = "")
  grDevices::png(
    loadingplot,
    width = 8,
    height = 8,
    units = "in",
    res = 300
  )
  graphics::plot(
    pc1L,pc2L,
    xlab=xlabeL,ylab=ylabeL,
    main = "Loading plot",
    pch = 19
  )
  grDevices::dev.off()
}
