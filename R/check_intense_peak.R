#' check the ppm error for the intense peaks
#'
#'
#' @export
check_intense_peak = function(RT_mass, RI_lib , RI_sample, n_peaks=1, acceptable_PPM_err = 10) {

  for(z in 1:length(RT_mass)){

    RT_mass[[z]]$intense_peak <- F
    RT_mass[[z]]$intense_peak2 <- NA
    RT_mass[[z]]$intense_peak3 <- NA


    z_peaks <- as.data.frame(RI_lib[names(RI_lib) == names(RT_mass)[z]])


    for(zi in 1: nrow(RT_mass[[z]])){

      zi_peaks <- as.data.frame(RI_sample[names(RI_sample) == RT_mass[[z]]$Feature_ID[zi]])


      #the first
      table_values_a <- sort(unique(z_peaks[,2]), decreasing = T)
      table_values_b <- sort(unique(zi_peaks[,2]), decreasing = T)


      for(ii in 1:min(n_peaks, length(table_values_a))){

        a = z_peaks[z_peaks[,2] == table_values_a[ii], ]
        if(length(table_values_b) < ii){
          break
        }
        b = zi_peaks[zi_peaks[,2] == table_values_b[ii], ]

        PPM_err <- calc_ppm_err(a,b)
        RT_mass[[z]][zi, 6+ii] <- any(PPM_err < acceptable_PPM_err)
        if(!RT_mass[[z]][zi, 6+ii]){
          break
        }
      }


    }

    RT_mass[[z]] = RT_mass[[z]][apply(RT_mass[[z]][, c("intense_peak", "intense_peak2", "intense_peak3")], 1, all, na.rm=TRUE), ]
  }

  #filter the RT_mass by deleting the empty data.frame;
  RT_mass = RT_mass[sapply(RT_mass, function(x) dim(x)[1]) > 0]

return(RT_mass)

}

