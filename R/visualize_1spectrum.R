#' Visualize a single MS/MS spectrum
#'
#' @param spectrum MS/MS spectrum as two-columns data.frame
#' @param name spctrum label
#' @param peak_number whether to show the peak number or not
#' @importFrom ggplot2 geom_hline geom_segment geom_point scale_color_manual scale_y_continuous theme ylab scale_linewidth_manual
#' @importFrom ggrepel geom_text_repel
#' @returns ggplot object
#' @export

visualize_1spectrum <- function(spectrum=NULL, name=NULL, peak_number=FALSE){

  mz <- RI <- NULL
  if(is.null(name) ){
    stop("Provide a name.\n")
  }

   set_theme(theme_science())
  
  plot_data <- spectrum
  colnames(plot_data)[1:2] <- c("mz", "RI")
  
  seq_lines <- seq(0, 100, by=10)
 
  p <- ggplot(plot_data, aes(x=mz, xend=mz, y=0, yend=RI)) +
    geom_hline(yintercept = seq_lines, lty=2, col="gray")+
    geom_segment(colour = "#2C3E50") +
    geom_point(aes(x=mz, y=RI)) +
    theme(axis.title.y = element_text(angle = 0, vjust = 0.50)) +
    ylab("RI")
  
  if(peak_number){
    p <- p + geom_text_repel(aes(x=mz, y=RI))
  }
  
  return(p)
  
}
