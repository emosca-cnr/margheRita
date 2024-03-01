#' Missing values imputation
#' 
#' the NA values are replaced with a fraction of the minimum value of that row (min_r), sampled in the interval \eqn{[min_r * a, min_r * b]}
#' 
#' @param mRList mRList object
#' @param seed set this seed with set.seed
#' @param a minimum increase factor
#' @param b maxmimum increase factor
#' #' @export

imputation <- function(mRList, seed=NULL, a=0.1, b=0.25) {
  
  
  #internal function
  imp_row <- function(row_min, n_na, a=0.1, b=0.25, n=100){
    
    ans <- seq(from = row_min * a, to = row_min * b, by = 1/n)
    ans <- sample(ans, n_na)
    return(ans)
    
  }
  
  if(!is.null(seed)){
    cat("Setting seed =", seed, ".\n")
    set.seed(seed)
  }
  
  na_val <- is.na(mRList$data)
  n <- min(100, max(rowSums(na_val))) #n is set possibly to the maximum number of na found in a row.
  
  #index of rows with NA
  idx_na_rows <- apply(na_val, 1, any)
  
  if(any(idx_na_rows)){
    
    #replacement
    mRList$data[idx_na_rows, ] <- apply(mRList$data[idx_na_rows, ], 1, function(x) replace(x = x, list = is.na(x), values = imp_row(min(x, na.rm = T), sum(is.na(x)), a=a, b=b, n=n)))
    
  }else{
    cat("No NAs: imputation not performed.\n")
  }
  return(mRList)
  
}

