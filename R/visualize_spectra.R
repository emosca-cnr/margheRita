#' Visualize the spectra of good candidates
#'
#' @param candidate_list annotation_list
#' @param x one or more metabolites to visualize
#' @export
#'
visualize_spectra <- function(annotation_list, x=NULL, out_dir="./"){
  
  for(z in 1:length(x)){ # cycle thought library metabolites
    
    candidates <- annotation_list$candidate_list[names(annotation_list$candidate_list) == x[z]][[1]]
    reference_peaks <- annotation_list$RI_reference[names(annotation_list$RI_reference) == x[z]][[1]]
    
    jpeg(file=paste0(out_dir, "/", x[z], ".jpg"), width = 200, height = 200, res=300, units = "mm")
    par(mfrow = c(round(nrow(candidates) / 2 + 1), round(nrow(candidates) / 2)))
    
    for(zi in 1:nrow(reference_peaks)){ #cycle through candidates
      
      candidate_peaks <- annotation_list$RI_data[names(annotation_list$RI_data) == candidates$Feature_ID[zi]][[1]]
      
      plot(reference_peaks[, 1] , reference_peaks[, 2] ,
           xlim=c(min(reference_peaks[,1], candidate_peaks[,1]), max(reference_peaks[,1], candidate_peaks[,1])), ylim = c(-100, 100),
           type = "h" , col= "red" ,
           xlab = names(candidate_peaks)[1] , ylab = "Relative Intensity",
           main= names(annotation_list$candidate_list[1]), font.main= 2 , cex.main= 0.9 )
      
      lines( candidate_peaks[,1], -candidate_peaks[,2], type = "h" , col="blue")
      
    }
    
    dev.off()
  }
  
  return()
}
