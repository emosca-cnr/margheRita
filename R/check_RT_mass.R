#' Merge candidates with appropriate retention time and PPM error
#'
#'
#' @export
#'
#' @param RT A list of library ID each contain a data frame of sample ID with ideal retention time. it is obtained from previous function.
#' @param mass A list of library ID each contain a data frame of sample ID with ideal PPM error. it is obtained from previous function.
#'
#' @return data.frame with candidates

check_RT_mass <- function(RT=NULL, mass=NULL){

  
  RT_mass <- intersect(names(RT), names(mass))
  RT_mass <- setNames(vector("list", length(RT_mass)), RT_mass)

  RT <- RT[match(names(RT_mass), names(RT))]
  mass <- mass[match(names(RT_mass), names(mass))]
  
  for (z in 1:length(RT_mass)) {

    RT_mass[[z]] <-  merge(RT[[z]], mass[[z]], by= "Feature_ID")
  }

  #filter the RT_mass by deleting the empty data.frame;
  RT_mass <- RT_mass[sapply(RT_mass, function(x) dim(x)[1]) > 0]
  
  RT_mass <- lapply(RT_mass, function(x) x[order(x$RT_err * x$ppm_error), ])

  return(RT_mass)
}
