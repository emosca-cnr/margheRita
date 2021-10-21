#' pareto scaling of the dataframe after missing value correction
# dataframe only with samples and intensity of each m/z and RT
#' @param df data.frame without the first three columns (MS DIAL, RT and m/z)
#' @importFrom utils write.csv
#' @export
#' @importFrom stats sd


pareto <- function(X, centering=TRUE){


  # Here we perform centering
  if(centering){
    X <- apply(X, 1, function(x) x - mean(x))
  }

  # Then we perform scaling on the mean-centered matrix
  X <- apply(X, 1, function(x) x/sqrt(sd(x)))

  return(X)
}




