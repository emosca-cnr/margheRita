#' check_mass
#'
#'
#' @export


check_mass <- function(feature_data, reference){

  feature_data$PEPMASS <- as.numeric(feature_data$PEPMASS)
  feature_data$RTINMINUTES <- as.numeric(feature_data$RTINMINUTES)

  mass = vector("list", nrow(reference))
  names(mass) = reference$`Metabolite name`

  for (j in 1:dim(reference)[1]) {

    ppm_error = abs((reference$`Average Mz`[j] - feature_data$PEPMASS)) / reference$`Average Mz`[j] * 100


    mass[[j]] = data.frame( ID = feature_data$ID , ppm_error= ppm_error,  stringsAsFactors = F)

    mass[[j]]$mass_flag = mass[[j]]$ppm_error< 15
    mass[[j]] = mass[[j]][mass[[j]]$mass_flag,]
    mass[[j]]$mass_status = "super"
    mass[[j]]$mass_status[mass[[j]]$ppm_error >= 5] = "acceptable"
    mass[[j]]$mass_status[mass[[j]]$ppm_error >= 10] = "suffer"

  }

  return(mass)
}
