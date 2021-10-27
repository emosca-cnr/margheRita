#' check_mass
#'
#'
#' @export


check_mass <- function(feature_data, reference){

   mass = vector("list", nrow(reference))
  names(mass) = reference$name

  for (j in 1:dim(reference)[1]) {

    ppm_error = abs((reference$mz[j] - feature_data$mz)) / reference$mz[j] * 100


    mass[[j]] = data.frame( Feature_ID = feature_data$Feature_ID , ppm_error= ppm_error,  stringsAsFactors = F)

    mass[[j]]$mass_flag = mass[[j]]$ppm_error< 15
    mass[[j]] = mass[[j]][mass[[j]]$mass_flag,]
    mass[[j]]$mass_status = "super"
    mass[[j]]$mass_status[mass[[j]]$ppm_error >= 5] = "acceptable"
    mass[[j]]$mass_status[mass[[j]]$ppm_error >= 10] = "suffer"

  }

  return(mass)
}
