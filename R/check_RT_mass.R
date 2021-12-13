#' Merging ideal candidate in term of Retention Time and PPM error
#'
#'
#' @export
#'
check_RT_mass = function( RT , mass, reference ){

  RT_mass = vector("list", nrow(reference))
  names(RT_mass) = reference$Name

  for (z in 1:dim(reference)[1]) {

    RT_mass[[z]] <-  merge(RT[[z]] , mass[[z]] , by= "Feature_ID")
  }

  #filter the RT_mass by deleting the empty data.frame;
  RT_mass = RT_mass[sapply(RT_mass, function(x) dim(x)[1]) > 0]

  return(RT_mass)
}
