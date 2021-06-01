#' missing values imputation
#' @param dataframe of input contains intensity of all features according to samples
#' @export

imputation <- function(m_list, seed=1234) {

  set.seed(seed) #hec

  for (i in 1:nrow(m_list$df)){
#
    min <- min(m_list$df[i, 4:ncol(m_list$df)],na.rm=T)
    min_val <- seq(from = (min * 0.1), to = (min * 0.25), by = (min * 0.01)) #if min*0.025 is greater than the row max? possible solution... to = min(max(m_list$df[i, 4:ncol(m_list$df)], na.rm=T), min * 0.25)

    if(any(is.na(m_list$df[i, 4:ncol(m_list$df)]))){
      idx_i <- is.na(m_list$df[i, 4:ncol(m_list$df)])
      m_list$df[i, idx_i] <- sample(min_val, size = sum(idx_i), replace = T)
    }
  }
  return(m_list)
}

