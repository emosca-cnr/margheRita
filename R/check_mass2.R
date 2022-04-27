#' selection of the mass for family of the metabolites. 
#' @export
#'
#'
check_mass2 <- function(feature_data = NULL, reference = NULL, unaccept_flag= 1000){
  
  mass = vector("list", nrow(reference))
  names(mass) = reference$ID
  
  
  for (j in 1:dim(reference)[1]) {
    
    ppm_error = abs((reference$mz[j] - feature_data$mz)) / reference$mz[j] * 1000000
    mass[[j]] = data.frame( Feature_ID = feature_data$Feature_ID , ppm_error= ppm_error,  stringsAsFactors = F)
    
    if(nrow(mass[[j]]) > 0){   #issue with empty data.frame
      
      mass[[j]]$mass_flag = mass[[j]]$ppm_error < unaccept_flag
      mass[[j]] = mass[[j]][mass[[j]]$mass_flag,]
    }
  }
  
  return(mass)
}