#' Calculate PPM error
#' 
#' PPM_error = abs(x - y) / x * 1000000
#'
#' @param a data.frame of 1 row and PPM in column 1
#' @param b data.frame of 1 row and PPM in column 1
#' @return PPM error

calc_ppm_err <- function(a=NULL, b=NULL){
  
  PPM_err <- matrix(0, nrow(a), nrow(b))
  
  for(i in 1:nrow(PPM_err)){
    for(j in 1:ncol(PPM_err)){
      PPM_err[i, j] = abs(a[i, 1] - b[j, 1]) / a[i, 1] * 1000000
      
    }
  }
  return(PPM_err)
}