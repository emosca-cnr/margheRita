#' check_mass
#'
#'
#' @export


check_mass <- function(feature_data, reference, unaccept_flag= unaccept_flag, accept_flag= accept_flag , suffer_flag= suffer_flag ){

   mass = vector("list", nrow(reference))
  names(mass) = reference$Name

  for (j in 1:dim(reference)[1]) {

    ppm_error = abs((reference$mz[j] - feature_data$mz)) / reference$mz[j] * 1000000

    mass[[j]] = data.frame( Feature_ID = feature_data$Feature_ID , ppm_error= ppm_error,  stringsAsFactors = F)

    if(nrow(mass[[j]]) > 0){   #issue with empty data.frame

    mass[[j]]$mass_status = "super"  #comes here because it can not add column to empty data.frame

    mass[[j]]$mass_flag = mass[[j]]$ppm_error < unaccept_flag
    mass[[j]] = mass[[j]][mass[[j]]$mass_flag,]


    mass[[j]]$mass_status[mass[[j]]$ppm_error >= accept_flag] = "acceptable"
    mass[[j]]$mass_status[mass[[j]]$ppm_error >= suffer_flag] = "suffer"

    }
  }

  #filter the mass by deleting the empty data.frame;
  #mass = mass[sapply(mass, function(x) dim(x)[1]) > 0]


  return(mass)
}
