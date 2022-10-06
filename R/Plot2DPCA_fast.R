#' Create 2D PCA score and loading plots
#' Create 2D PCA score and loading plots choosing the principal components.
#' @param mRList mRList object
#' @param pcx Specify the principal component on the x-axis
#' @param pcy Specify the principal component on the y-axis
#' @param dirout output directory
#' @param col_by (default class)
#' @param include_QC (default TRUE)
#' @export
#' @importFrom pcaMethods scores
#' @importFrom graphics plot legend
#' @importFrom grDevices png rainbow dev.off

Plot2DPCA_fast <- function(mRList, pcx=1, pcy=2, dirout = "./", col_by="class", include_QC=TRUE){
  
  dir.create(dirout)
  
  xlabel = paste("PC", pcx, "(", mRList$pca@R2[pcx], "%)")
  ylabel = paste("PC", pcy, "(", mRList$pca@R2[pcy], "%)")
  xy <- scores(mRList$pca)
  xy <- xy[, c(pcx, pcy)]
  
  
  scoreplot = paste(dirout, "/Scoreplot.png", sep = "")
  
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
  
  grDevices::png(scoreplot, width = 200, height = 200, units = "mm", res = 300)
  
  graphics::plot(
    xy[, 1], xy[, 2],
    xlab=xlabel, ylab=ylabel,
    main = "Score plot",
    col = col_pal[as.numeric(col_factor)],
    pch = 19
  )
  legend("bottomright", legend = levels(col_factor), col = col_pal, pch=16, cex=0.5)
  
  dev.off()

 }
