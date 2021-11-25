#'calculating and selecting proper PPM error and RI for mass mass
#'
#'
#' @export
#'
check_mass_mass = function(RT_mass, RI_lib , RI_sample) {


  for(z in 1: length(RT_mass)){

    RT_mass[[z]]$peaks_found_ppm_RI <- 0
    RT_mass[[z]]$peaks_found_flag <- "red"


    z_peaks <- as.data.frame(RI_lib[names(RI_lib) == names(RT_mass)[1]])
    #if (length(RI_lib[[z]]) <= 2) { z_peaks = t(z_peaks)}

    #compare ppm error of peaks between z and every candidate z_i
    for(zi in 1:nrow(RT_mass[[z]])){
      zi_peaks <- as.data.frame(RI_sample[names(RI_sample) == RT_mass[[z]]$Feature_ID[zi]])
      #if (length(RI_lib[[z]]) <= 2) { zi_peaks = t(zi_peaks)}

      #calculate ppm error between the peaks of z and the peaks of z_i
      pmm_error_matrix <- matrix(0, nrow(z_peaks), nrow(zi_peaks))
      pmm_error_matrix_flags <- pmm_error_matrix

      RI = pmm_error_matrix
      RI_flags = pmm_error_matrix


      for(i in 1:nrow(pmm_error_matrix)){
        for(j in 1:ncol(pmm_error_matrix)){
          pmm_error_matrix[i, j] = abs(z_peaks[i, 1] - zi_peaks[j, 1]) / z_peaks[i, 1] * 100

          RI[i, j] = abs(z_peaks[i, 2] - zi_peaks[j, 2])
        }
      }

      #count the peaks that are found
      pmm_error_matrix_flags[pmm_error_matrix <10] <- 1
      RI_flags[RI < 30] <- 2

      flags_matrix <- pmm_error_matrix_flags + RI_flags

      if(nrow(flags_matrix) > ncol(flags_matrix)) {
        RT_mass[[z]]$peaks_found_ppm_RI[zi] <- sum(sign(colSums(flags_matrix==3)))
      }

        if(nrow(flags_matrix) < ncol(flags_matrix)) {
          RT_mass[[z]]$peaks_found_ppm_RI[zi] <- sum(sign(rowSums(flags_matrix==3)))
        }
  }

    #define the flags
    RT_mass[[z]]$peaks_found_flag[RT_mass[[z]]$peaks_found_ppm_RI > 1] <- "yellow"
    RT_mass[[z]]$peaks_found_flag[RT_mass[[z]]$peaks_found_ppm_RI > 2] <- "green"
}

    #filtering RT_mass by deleting the candidates with bad ppm error and RI
    for (k in 1:length(RT_mass)) {
      RT_mass[[k]] = RT_mass[[k]][RT_mass[[k]]$peaks_found_flag != "red", ]
    }


  return(RT_mass)
}
