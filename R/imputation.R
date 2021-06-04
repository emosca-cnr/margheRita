#' missing values imputation
#' @param dataframe of input contains intensity of all features according to samples
#' @export

imputation <- function(m_list, seed=NULL, a=0.1, b=0.25, n=100) {


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

  na_val <- is.na(m_list$data)
  n <- min(n, max(rowSums(na_val))) #n is set possibly to the maximum number of na found in a row.

  #index of rows with NA
  idx_na_rows <- apply(na_val, 1, any)
  if(any(idx_na_rows)){

    #replacement
    m_list$data[idx_na_rows, ] <- apply(m_list$data[idx_na_rows, ], 1, function(x) replace(x = x, list = is.na(x), values = imp_row(min(x, na.rm = T), sum(is.na(x)), a=a, b=b, n=n)))

  }else{
    cat("Nothing to do.\n")
  }
  return(m_list)

}

