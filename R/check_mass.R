#' calculating PPM error
#'
#'`check_mass` it is making a list of all library metabolites and assigning proper sample ID to each based on desired PPM error.
#'
#' @export
#'
#' @param reference A list of library contain mz(mass-to-charge ratio) with specific ID.
#' @param feature_data A list of sample data contain mz(mass-to-charge ratio) with specific ID.
#' @param accept_flag A number with default value of 5. PPM errors < accept_flag will be tagged as "super", while those > accept_flag and < suffer_flag will be tagged as "acceptable"
#' @param suffer_flag A number with default value of 10. PM errors above this value and < unaccept_flag will be tagged as "suffer"
#' @param unaccept_flag A number with default value of 15. The maximum PPM error must be less than this value. and those above this number will be eliminated.

#' @return mass, A list of library ID each contain a data frame of sample ID with a range of PPM error less than unacceptable flag
#' @examples


check_mass <- function(reference=NULL, feature_data=NULL, unaccept_flag=15, accept_flag=5, suffer_flag=10){

   mass = vector("list", nrow(reference))
  #names(mass) = reference$Name
   names(mass) = reference$ID

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

  return(mass)
}
