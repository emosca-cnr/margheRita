#'visualize the spectra of good candidates
#'
#'
#' @export
#'
visualize_spectra2 <- function(x=NULL, RI_sample=NULL, RI_lib=NULL){

  for(a in 1: length(x)){ # cycle thought library metabolites
    lib_plot <- as.data.frame(RI_lib[names(RI_lib) == x[[a]]$ID_peaks[a]])

    jpeg(file=paste0(dir_out, "/", ".jpg"))

if(nrow(x[[a]])>1){
    par(mfrow = c(round(nrow(x[[a]])/2), round(nrow(x[[a]])/2)+1))
}else{par(mfrow = c(round(nrow(x[[a]])/2)+1, round(nrow(x[[a]])/2)+1))}

    for(b in 1:nrow(x[[a]])){ #cycle through candidates
      sample_plot <- as.data.frame(RI_sample[names(RI_sample) == x[[a]]$Feature_ID[b]])

if(length(lib_plot)>0){
      plot(lib_plot[,1] , lib_plot[,2] ,
           xlim=c(min(lib_plot[ ,1], sample_plot[ ,1]), max(lib_plot[ ,1], sample_plot[ ,1])),
           ylim = c(-100, 100),
           type = "h" , col= "red" ,
           xlab = names(sample_plot)[1] , ylab = "Relative Intensity",
           main= unique(x[[a]]$Name), font.main= 2 , cex.main= 0.9 )

      lines( sample_plot[ ,1], -sample_plot[ ,2], type = "h" , col="blue")
}
    }
    dev.off()
  }
return()

}
