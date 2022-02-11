#' equation to calculate the PPM error
#'
#'
#' @export


calc_ppm_err <- function(a=NULL, b=NULL){
  
  PPM_err <- matrix(0, nrow(a), nrow(b))
  
  for(i in 1:nrow(PPM_err)){
    for(j in 1:ncol(PPM_err)){
      PPM_err[i, j] = abs(a[i, 1] - b[j, 1]) / a[i, 1] * 1000000
      
    }
  }
  return(PPM_err)
}