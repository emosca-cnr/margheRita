#' Visualize feature-metabolite MS/MS spectra
#'
#' @param spectra_top MS/MS spectrum as two-columns data.frame
#' @param spectra_bottom MS/MS spectrum as two-columns data.frame
#' @param name_top name for top spectrum
#' @param name_bottom name for bottom spectrum
#' @param type mirrored or overlapped visualization
#' @param peak_number whether to show the peak number or not
#' @param match_matrix match matrix calculated by metabolite_identi
#' @importFrom ggplot2 geom_hline geom_segment geom_point scale_color_manual scale_y_continuous theme ylab scale_linewidth_manual
#' @importFrom ggrepel geom_text_repel
#' @returns ggplot object
#' @export

visualize_feature_metabolite_spectra <- function(spectra_top=NULL, spectra_bottom=NULL, name_top=NULL, name_bottom=NULL, type=c("mirrored", "overlapped"), peak_number=FALSE, match_matrix=NULL){

  mz <- RI <- spectrum <- label <- matches <- NULL
  if(is.null(name_top) | is.null(name_bottom)){
    stop("Provide names for top and bottom.\n")
  }
  type <- match.arg(type)

  peak_matches <- rep(FALSE, nrow(spectra_top)+nrow(spectra_bottom))
  if(!is.null(match_matrix)){
    lib_n <- which(apply(match_matrix, 1, function(x) any(x==3)))
    feat_n <- which(apply(match_matrix, 2, function(x) any(x==3)))
    if(length(lib_n)>0 & length(feat_n)>0){
      peak_matches[lib_n] <- TRUE
      peak_matches[nrow(spectra_top)+feat_n] <- TRUE
    }
  }
  
  set_theme(theme_science())
  
  plot_data <- data.frame(rbind(spectra_top, spectra_bottom), spectrum=c(rep(name_top, nrow(spectra_top)), rep(name_bottom, nrow(spectra_bottom))), label=c(1:nrow(spectra_top), 1:nrow(spectra_bottom)), matches=peak_matches)
  colnames(plot_data)[1:2] <- c("mz", "RI")
  
  seq_lines <- seq(0, 100, by=10)
  if(type=="mirrored"){
    plot_data$RI[plot_data$spectrum == name_bottom] <- -plot_data$RI[plot_data$spectrum == name_bottom]
    seq_lines <- seq(-100, 100, by=10)
  }
  
  p <- ggplot(plot_data, aes(x=mz, xend=mz, y=0, yend=RI, col=spectrum, label=label, linewidth = matches)) +
    geom_hline(yintercept = seq_lines, lty=2, col="gray")+
    geom_segment() +
    geom_point(aes(x=mz, y=RI)) +
    scale_linewidth_manual(values = c(0.5, 1.5)) +
    scale_color_manual(values=c("#E05A47", "#2C3E50")) +
    scale_y_continuous(labels = function(x) gsub("-", "", x)) +
    theme(axis.title.y = element_text(angle = 0, vjust = 0.50)) +
    ylab("RI")
  
  if(peak_number){
    p <- p + geom_text_repel(aes(x=mz, y=RI), show.legend = FALSE)
  }
  
  return(p)
  
}
