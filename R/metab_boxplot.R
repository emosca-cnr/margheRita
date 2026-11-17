#' Boxplot of feature abundance
#'
#' @param mRList mRList object
#' @param feature feature id
#' @param col_by define how to color the boxplot, default is class
#' @param col_pal optional palette. It must match the number of levels of the factor indicated by col_by
#' @param show_stat whether to show the p-values from t test, adjusted by FDR
#' @importFrom ggplot2 set_theme ggplot geom_boxplot stat_boxplot scale_fill_manual theme xlab ylab
#' @importFrom pals alphabet
#' @importFrom ggpubr stat_compare_means geom_pwc
#' @importFrom rlang .data
#' @returns ggplot object
#' @export

metab_boxplot <- function(mRList=NULL, feature=NULL, col_by="class", col_pal=NULL, show_stat=FALSE){
  
  y <- NULL
  X_ann <- mRList$sample_ann
  col_factor <- as.factor(X_ann[, col_by])
  col_factor_lev <- levels(col_factor)
  if(is.null(col_pal) | length(col_pal) < length(col_factor_lev)){
    cat("Using default palette\n")
    col_pal <- as.character(alphabet(length(col_factor_lev)))
  }
  
  
  data <- mRList$data[rownames(mRList$data) %in% feature, , drop=FALSE] #select the metabolites (feature), according to list from users
  
  if(nrow(data)<1){
    message("The feature was not found in the dataset.\n")
  }
  
  plot_data <- data.frame(id=colnames(mRList$data), y=as.numeric(mRList$data[rownames(mRList$data) == feature, ]))
  plot_data <- merge(plot_data, mRList$sample_ann, by.x="id", by.y=0)
  
  
  set_theme(theme_science())
  
  p <- ggplot(plot_data, aes(x=.data[[col_by]], y=y, fill=.data[[col_by]])) +
    geom_boxplot(outliers = F, lwd=.2) +
    stat_boxplot(geom='errorbar', linetype=1, width=0.5, lwd=.2) +
    scale_fill_manual(values = col_pal) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(y="abundance", x="", title = feature) +
    theme(plot.title = element_text(hjust = 0.5))
  
  if(show_stat){
    p <- p + geom_pwc(aes(group = .data[[col_by]]), method = "t_test", label = "p.format", label.size = 3, p.adjust.method = "fdr")
  }
    
  return(p)
  
}
