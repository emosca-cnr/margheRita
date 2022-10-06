#' Pareto scaling
#' Perform data scaling with Pareto approach.
#' @param X is dataframe that contains only data, without the columns that contain ID metabolite, m/z, RT, MS and MS/MS information
#' @param centering whether to center the data or not. Defaults to TRUE.
#' @importFrom utils write.csv
#' @export
#' @importFrom stats sd


pareto <- function(X, centering=TRUE){


  # Here we perform centering
  if(centering){
    X <- t(apply(X, 1, function(x) x - mean(x))) #ettore: t is necessary to keep a features-by-samples 
  }

  # Then we perform scaling on the mean-centered matrix
  X <- t(apply(X, 1, function(x) x/sqrt(sd(x)))) #ettore: t is necessary to keep a features-by-samples 

  
  
  return(X)
}




