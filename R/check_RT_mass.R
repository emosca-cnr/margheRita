#' Merge candidates with appropriate retention time and PPM error
#'
#'
#' @export
#'
#' @param RT A list of library ID each contain a data frame of sample ID with ideal retention time. it is obtained from previous function.
#' @param mass A list of library ID each contain a data frame of sample ID with ideal PPM error. it is obtained from previous function.
#' @param reference A list of library contain retention time and mz (mass-to-charge ratio) with specific ID.
#'
#' @return data.frame with candidates

check_RT_mass <- function(RT=NULL, mass=NULL, reference=NULL ){

  RT_mass = vector("list", nrow(reference))
  #names(RT_mass) = reference$Name
  names(RT_mass) = reference$ID

  for (z in 1:dim(reference)[1]) {

    RT_mass[[z]] <-  merge(RT[[z]] , mass[[z]] , by= "Feature_ID")
  }

  #filter the RT_mass by deleting the empty data.frame;
  RT_mass <- RT_mass[sapply(RT_mass, function(x) dim(x)[1]) > 0]
  
  RT_mass <- lapply(RT_mass, function(x) x[order(x$RT_err*x$ppm_error), ])

  return(RT_mass)
}
