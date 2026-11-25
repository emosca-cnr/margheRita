#' Relative log abundance plot
#' @param mRList mRList object
#' @param col_by (default class)
#' @param include_QC (default FALSE)
#' @param xlabels whether the show x axis labels
#' @param col_pal optional palette. It must match the number of levels of the factor indicated by #' @return  ggplot2 object
#' @importFrom ggplot2 set_theme ggplot aes geom_boxplot theme element_text labs geom_hline stat_boxplot
#' @importFrom reshape2 melt
#' @importFrom rlang .data
#' @returns ggplot object
#' @export 

plot_RLA <- function(mRList=NULL, col_by="class", col_pal=NULL, include_QC=FALSE, xlabels=FALSE){
  
  Var2 <- value <- id <- NULL
  
  if(include_QC){
    sample_ann <- rbind(mRList$sample_ann, mRList$QC_ann)
  }else{
    sample_ann <- mRList$sample_ann
  }
  
  stopifnot(length(mRList$RLA)>0)
  stopifnot(ncol(mRList$RLA) == nrow(sample_ann))
  
  lgd <- FALSE
  if(!is.null(col_by)){
    stopifnot(col_by %in% colnames(sample_ann))
    lgd <- TRUE
    if(is.null(col_pal)){
      col_pal <- as.character(pals::alphabet(length(levels(factor(sample_ann[, col_by])))))
    }
  }else{
    col_pal <- "gray"
    col_by <- "fill"
  }
  
  plot_data <- melt(mRList$RLA)
  plot_data <- merge(plot_data, sample_ann, by.x="Var2", by.y=0)
  if(col_by == "fill"){
    plot_data$fill <- ""
  }
  
  set_theme(theme_science())
  
  p <- ggplot(plot_data, aes(x=Var2, y=value, fill=.data[[col_by]], group=id)) +
    geom_boxplot(outliers = F, lwd=.2) +
    stat_boxplot(geom='errorbar', linetype=1, width=0.5, lwd=.2) +
    scale_fill_manual(values = col_pal) +
    geom_hline(yintercept=0, linetype=2, lwd=.2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.y = element_text(angle = 0, vjust = 0.5)) +
    labs(x="", y="RLA")
  if(!lgd){
    p <- p + theme(legend.position = "none")
  }
  if(!xlabels){
    p <- p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }
  
  return(p)
  
}
