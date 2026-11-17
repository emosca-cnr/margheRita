#' Calculate the PPM error of precursors
#' Calculate the PPM error of precursors and assign a flag according to its value.
#'
#' @param reference data frame of reference m/z values
#' @param feature_data data frame of m/z values observed in the samples 
#' @param accept_flag A number with default value of 5. PPM errors < accept_flag will be tagged as "super", while those > accept_flag and < suffer_flag will be tagged as "acceptable"
#' @param suffer_flag A number with default value of 10. PM errors above this value and < unaccept_flag will be tagged as "suffer"
#' @param unaccept_flag A number with default value of 15. The maximum PPM error must be less than this value. and those above this number will be eliminated.
#' @param filter whether to remove unacceptable associations

#' @return A list of where each element is named with a reference library ID and contains a data frame of matched features


check_mass <- function(reference=NULL, feature_data=NULL, unaccept_flag=20, accept_flag=5, suffer_flag=10, filter=TRUE){
  
  mass = vector("list", nrow(reference))
  names(mass) = reference$ID
  
  for (j in 1:dim(reference)[1]) {
    
    ppm_error = abs((reference$mz[j] - feature_data$mz)) / reference$mz[j] * 1000000
    
    mass[[j]] = data.frame(Feature_ID = feature_data$Feature_ID, mz=feature_data$mz, ppm_error= ppm_error,  stringsAsFactors = F)
    
    if(nrow(mass[[j]]) > 0){   #issue with empty data.frame
      
      mass[[j]]$mass_flag = mass[[j]]$ppm_error <= unaccept_flag
      
      mass[[j]]$mass_status = "super"  #comes here because it can not add column to empty data.frame
      mass[[j]]$mass_status[mass[[j]]$ppm_error > accept_flag] = "acceptable"
      mass[[j]]$mass_status[mass[[j]]$ppm_error > suffer_flag] = "suffer"
      mass[[j]]$mass_status[mass[[j]]$ppm_error > unaccept_flag] = "unacceptable"
      
      if(filter){
        mass[[j]] = mass[[j]][mass[[j]]$mass_flag, ]
      }
      
    }
  }
  
  #remove library elements without any match
  mass <- mass[unlist(lapply(mass, function(x) nrow(x)>0))]
  
  return(mass)
}
