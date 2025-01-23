#' Plot score and loading plot
#' @description Plot 2D PCA score and loading plots choosing the principal components.
#' @param mRList mRList object
#' @param pcx Specify the principal component on the x-axis
#' @param pcy Specify the principal component on the y-axis
#' @param dirout output directory
#' @param col_by (default class)
#' @param include_QC (default TRUE)
#' @param name_samples whether the names of samples should be added to each point
#' @param name_samples_cex if name_samples is TRUE, specify here the size of the sample names 
#' @param col_pal optional palette. It must match the number of levels of the factor indicated by col_by
#' 
#' @export
#' @importFrom pcaMethods scores
#' @importFrom grDevices png dev.off
#' @importFrom graphics layout plot.new text
#' @importFrom pals alphabet
#' @importFrom plotrix thigmophobe.labels

Plot2DPCA <- function(mRList=NULL, pcx=1, pcy=2, dirout = NULL, col_pal=NULL, col_by="class", include_QC=FALSE, name_samples = FALSE, name_samples_cex = 0.8){

  if (!is.null(dirout)) {
    dir.create(dirout, showWarnings = F)
  }else{
    dirout <- getwd()
  }

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
  col_factor_lev <- levels(col_factor)
  if(is.null(col_pal) | length(col_pal) < length(col_factor_lev)){
    cat("Using default palette\n")
    col_pal <- alphabet(length(col_factor_lev))
  }

  png(file.path(dirout, "scores.png"), width = 200, height = 160, units = "mm", res = 300)

  layout(matrix(c(1:2), nrow = 1), widths = c(.85, .15))
  
  par(mar=c(3, 3, .5, .5))
  par(mgp=c(1.6, .5, 0))
  
  plot(
    xy[, 1], xy[, 2],
    xlab=xlabel, ylab=ylabel,
    main = "",
    col = col_pal[as.numeric(col_factor)],
    pch = 19
  )
  if (name_samples) {
    #text(x = xy[, 1], y = xy[, 2], labels = rownames(xy), cex = name_samples_cex, pos=3, xpd=TRUE)
    thigmophobe.labels(x = xy[, 1], y = xy[, 2], labels = rownames(xy), cex = name_samples_cex)
  }
  
  
  par(mar=c(.1, .1, .1, .1))
  plot.new()
  legend("center", legend = levels(col_factor), col = col_pal, pch=16, cex=.6)

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
  if (name_samples) {
    #text(xyL[,1], xyL[,2], labels=rownames(xyL), cex=name_samples_cex, pos=3, xpd=TRUE)
    thigmophobe.labels(xyL[,1], xyL[,2], labels=rownames(xyL), cex=name_samples_cex)
  }
  

  dev.off()
}
