#' peak matching
#'
#' @export
#'
peak_matching = function(feature_data=NULL, reference=NULL, RT_mass=NULL, RI_lib=NULL, RI_sample=NULL,  lib_peaks_data=NULL, mode=c("POS", "NEG"), ppm_err=10, intensity=30) {

  #a novel data to not touch the input
  ans <- vector("list", length(RT_mass))
  names(ans) <- names(RT_mass)

  #ensure the same order
  lib_peaks_data <- lib_peaks_data[match(names(RI_lib), lib_peaks_data$ID), ]

  for(z in 1:length(RT_mass)){
    RT_mass[[z]]$peaks_found_ppm_RI <- 0
    RT_mass[[z]]$ID <- names(RT_mass)[z]

    CAS_z <- reference$CAS[reference$ID == names(RT_mass)[z]]
    z_peaks_all <- RI_lib[lib_peaks_data$CAS == CAS_z & lib_peaks_data$Collision_energy == mode ]

    if(length(z_peaks_all)>0){

      for(zzzz in 1:length(z_peaks_all)){
        z_peaks <- as.data.frame(z_peaks_all[[zzzz]])

        #RT_mass[[z]]$ID_peaks <- names(z_peaks_all)[zzzz]


        for(zi in 1:nrow(RT_mass[[z]])){
          zi_peaks <- as.data.frame(RI_sample[names(RI_sample) == RT_mass[[z]]$Feature_ID[zi]])

          pmm_error_matrix <- matrix(0, nrow(z_peaks), nrow(zi_peaks))
          pmm_error_matrix_flags <- pmm_error_matrix

          RI <- pmm_error_matrix
          RI_flags <- pmm_error_matrix

          for(i in 1:nrow(pmm_error_matrix)){
            for(j in 1:ncol(pmm_error_matrix)){
              pmm_error_matrix[i, j] = abs(z_peaks[i, 1] - zi_peaks[j, 1]) / z_peaks[i, 1] * 1000000

              RI[i, j] = abs(z_peaks[i, 2] - zi_peaks[j, 2])
            }
          }

          #count the peaks that are found
          pmm_error_matrix_flags[pmm_error_matrix < ppm_err] <- 1
          RI_flags[RI < intensity] <- 2

          flags_matrix <- pmm_error_matrix_flags + RI_flags


          RT_mass[[z]]$peaks_found_ppm_RI[zi] <- sum(sign(rowSums(flags_matrix==3)))



          if(is.null(ans[[z]])){
            ans[[z]] <- data.frame(Feature_ID=RT_mass[[z]]$Feature_ID[zi], ID_peaks=names(z_peaks_all)[zzzz], stringsAsFactors = F)
          }else{
            ans[[z]] <- rbind(ans[[z]], data.frame(Feature_ID=RT_mass[[z]]$Feature_ID[zi], ID_peaks=names(z_peaks_all)[zzzz], stringsAsFactors = F))

          }

        }
      }
    }
    if(!is.null(ans[[z]])){
      #keep only peaks with TRUE flags
      ans[[z]] <- merge(RT_mass[[z]], ans[[z]], by="Feature_ID", sort=F,no.dups = TRUE, all.y=FALSE)
      ans[[z]] <- merge(ans[[z]], reference, by="ID")

    }

  }

  #filter the ans by deleting the empty data.frame;
  ans <- ans[sapply(ans, function(x) !is.null(x))]

  #filter the ans by deleting the candidates without match
  for (k in 1:length(ans)) {
    ans[[k]] <- ans[[k]][ans[[k]]$peaks_found_ppm_RI != 0, ] }

  #filter the ans by deleting the empty data.frame
  ans <- ans[sapply(ans, function(x) dim(x)[1]) > 0]

  return(ans)

}
