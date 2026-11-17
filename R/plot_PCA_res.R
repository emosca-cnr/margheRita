#' Generate all-pairs PC scatters and R2 explained by PCs
#' @param mRList mRList object
#' @param col_by (default class)
#' @param include_QC (default FALSE)
#' @param col_pal optional palette. It must match the number of levels of the factor indicated by col_by
#' 
#' @importFrom ggplot2 geom_point scale_color_manual theme ylab xlab geom_bar
#' @importFrom ggrepel geom_text_repel
#' @importFrom GGally ggpairs
#' @importFrom pcaMethods scores
#' @importFrom pals alphabet
#' @importFrom rlang .data
#' @returns List of two ggplot objects
#' @export

plot_PCA_res <- function(mRList=NULL, col_pal=NULL, col_by="class", include_QC=FALSE){

  PC <- R2 <- NULL
  
  #include QC
  X_ann <- mRList$sample_ann
  if(include_QC){
    cat("Including QC\n")
    X_ann<- rbind(X_ann, mRList$QC_ann)
  }
  
  #ensure the same order between annotation and data
  X_ann <- X_ann[match(rownames(mRList$pca@scores), X_ann[, "id"]), ]
  
  col_factor <- as.factor(X_ann[, col_by])
  col_factor_lev <- levels(col_factor)
  if(is.null(col_pal) | length(col_pal) < length(col_factor_lev)){
    cat("Using default palette\n")
    col_pal <- as.character(alphabet(length(col_factor_lev)))
  }
  
  set_theme(theme_science())
  p_list <- setNames(list("vector", length=2), c("scree", "pairs"))
  
  plot_data <- cbind(mRList$pca@scores, X_ann)
  p_list$pairs <-  ggpairs(plot_data, columns = 1:ncol(mRList$pca@scores), aes(colour=.data[[col_by]]))
  
  plot_data <- data.frame(PC=names(mRList$pca@R2), R2=as.numeric(mRList$pca@R2))
  p_list$scree <- ggplot(plot_data, aes(x=PC, y=R2)) +
    geom_bar(stat = "identity", fill="#2C3E50", color="#e9ecef",) +
    xlab("") +
    theme(axis.title.y = element_text(angle=0, vjust=0.5))
  

  return(p_list)

}
