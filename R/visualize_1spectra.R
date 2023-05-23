#' #' Visualize the spectra of good candidates
#' #'
#' #' @param spectra_top MS/MS spectra as two-columns data.frame
#' #' @param spectra_bottom MS/MS spectra as two-columns data.frame
#' #' @param name_top name for top spectra
#' #' @param name_bottomname for bottom spectra
#' #' @param peak_id peak id
#' #' @param out_file output file
#' #' @export
#' #' @importFrom graphics par lines legend
#' #' @importFrom grDevices jpeg adjustcolor
#' #' @importFrom plotrix thigmophobe.labels
#' #' @importFrom Hmisc minor.tick
#' #' @param type mirrored or overlapped visualization
#' #' @param peak_number whether to show the peak number or not
#' 
#' visualize_1spectra <- function(spectra_top=NULL, spectra_bottom=NULL, name_top=NULL, name_bottom=NULL, out_file=NULL, type=c("mirrored", "overlapped"), peak_number=FALSE){
#'   
#'   if(is.null(name_top) | is.null(name_bottom)){
#'     stop("Provide names for top and bottom.\n")
#'   }
#'   type <- match.arg(type)
#'   
#'   if(!is.null(out_file)){
#'     png(filename = out_file, width = 200, height = 200, res=300, units = "mm")
#'   }
#'   
#'   if(type=="mirrored"){
#'     ylim <- c(-100, 100)
#'     hlines <- seq(-100, 100, 10)
#'   }else{
#'     ylim <- c(0, 100)
#'     hlines <- seq(0, 100, 10)
#'   }
#'   
#'   #empty plot
#'   plot(0, pch="", xlim=c(min(spectra_top[, 1], spectra_bottom[, 1]), max(spectra_top[, 1], spectra_bottom[, 1])), ylim = ylim, xlab = "m/z", ylab = "Relative Intensity")
#'   abline(h=hlines, lty=2, col="gray")
#'   minor.tick(ny = 2, nx=1)
#'   
#'   #library
#'   points(spectra_top[, 1] , spectra_top[, 2], type = "h" , col= adjustcolor("red", 0.6), lwd=2)
#'   points(spectra_top[, 1], spectra_top[,2], col=adjustcolor("red", 0.8), pch=16)
#'   
#'   if(peak_number){
#'     plotrix::thigmophobe.labels(spectra_top[, 1], spectra_top[,2], 1:nrow(spectra_top))
#'   }
#'   
#'   ### feature
#'   if(type=="mirrored"){
#'     
#'     lines(spectra_bottom[, 1], -spectra_bottom[,2], type = "h" , col=adjustcolor("blue", 0.6), lwd=2)
#'     points(spectra_bottom[, 1], -spectra_bottom[,2], col=adjustcolor("blue", 0.6), pch=16)
#'     if(peak_number){
#'       plotrix::thigmophobe.labels(spectra_bottom[, 1], -spectra_bottom[, 2], 1:nrow(spectra_bottom))
#'     }
#'     
#'   }else{
#'     
#'     lines(spectra_bottom[, 1], spectra_bottom[,2], type = "h" , col=adjustcolor("blue", 0.6), lwd=2)
#'     points(spectra_bottom[, 1], spectra_bottom[,2], col=adjustcolor("blue", 0.6), pch=16)
#'     if(peak_number){
#'       plotrix::thigmophobe.labels(spectra_bottom[, 1], spectra_bottom[, 2], 1:nrow(spectra_bottom))
#'     }
#'     
#'   }
#'   
#'   legend("bottomright", legend=c(name_top, name_bottom), lty=1, col=c("red", "blue"), cex=0.7, bg="transparent") 
#'   
#'   if(!is.null(out_file)){
#'     dev.off()
#'   }
#'   
#' }
