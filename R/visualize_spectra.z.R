#'visualize the spectra of good candidates
#'
#'
#' @export
#'
visualize_spectra = function(annotation, RI_lib, RI_sample ) {

  for(z in 1: length(annotation)){ # cycle thought library metabolites
    z_peaks <- as.data.frame(RI_lib[names(RI_lib) == names(annotation)[z]])


    par(mfrow = c(round(nrow(annotation[[z]])/2),round(nrow(annotation[[z]])/2)))

    for(zi in 1:nrow(annotation[[z]])){ #cycle through candidates
      zi_peaks <- as.data.frame(RI_sample[names(RI_sample) == annotation[[z]]$Feature_ID[zi]])



      #jpeg(file= "")
      spectra = plot(z_peaks[,1] , z_peaks[,2] ,
           xlim=c(min(z_peaks[,1], zi_peaks[,1]), max(z_peaks[,1], zi_peaks[,1])), ylim = c(-100, 100),
           type = "h" , col= "red" ,
           xlab = names(zi_peaks)[1] , ylab = "Relative Intensity",
           main= names(annotation[1]), font.main= 2 , cex.main= 0.9 )


      lines( zi_peaks[,1], -zi_peaks[,2], type = "h" , col="blue")

      #dev.off()

    }
  }

  return(list(spectra))
}
