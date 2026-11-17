#' Plot scores and loadings
#' @description Plot 2D PCA scores and loadings choosing which pairs of PCs.
#' @param mRList mRList object
#' @param pcx Specify the principal component on the x-axis
#' @param pcy Specify the principal component on the y-axis
#' @param col_by (default class)
#' @param include_QC (default FALSE)
#' @param name_samples whether the names of samples should be added to each point
#' @param col_pal optional palette. It must match the number of levels of the factor indicated by col_by
#' @importFrom ggplot2 geom_point scale_color_manual theme ylab xlab
#' @importFrom ggrepel geom_text_repel
#' @importFrom pcaMethods scores
#' @importFrom pals alphabet
#' @importFrom rlang .data
#' @returns List of two ggplot objects
#' @export

Plot2DPCA <- function(mRList=NULL, pcx=1, pcy=2, col_pal=NULL, col_by="class", include_QC=FALSE, name_samples = FALSE){

  id <- NULL
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
  X_ann <- X_ann[match(rownames(xy), X_ann[, "id"]), ]

  col_factor <- as.factor(X_ann[, col_by])
  col_factor_lev <- levels(col_factor)
  if(is.null(col_pal) | length(col_pal) < length(col_factor_lev)){
    cat("Using default palette\n")
    col_pal <- as.character(alphabet(length(col_factor_lev)))
  }

  xlabeL = paste("Loading", pcx)
  ylabeL = paste("Loading", pcy)
  xyL <- as.data.frame(mRList$pca@loadings)
  xyL <- xyL[, c(pcx, pcy)]

  
  set_theme(theme_science())
  
  plot_data <- cbind(xy, X_ann)
  
  p_list <- setNames(list("vector", length=2), c("scores", "loadings"))
  
  p_list[[1]] <- ggplot(plot_data, aes(x=plot_data[, 1], y=plot_data[, 2], col=.data[[col_by]], label=id)) +
    geom_point() +
    scale_color_manual(values=as.character(col_pal)) + 
    xlab(xlabel) +
    ylab(ylabel)
  
  if (name_samples) {
    p_list[[1]] <- p_list[[1]] + geom_text_repel()
  }
  
  p_list[[2]] <- ggplot(xyL, aes(x=xyL[, 1], y=xyL[, 2])) +
    geom_point() +
    xlab(xlabeL) +
    ylab(ylabeL)
  
  return(p_list)

}
