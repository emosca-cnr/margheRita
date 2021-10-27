#' check_RT
#'
#'
#' @export


check_RT <- function(feature_data, reference){

  
  res <- vector("list", nrow(reference))
  names(res) <- reference$Name

  for (k in 1:dim(reference)[1]){

    idx_ok <- which((feature_data$rt > (reference$rt[k] - 1)) & (feature_data$rt < (reference$rt[k] +1)))

    res[[k]] <- data.frame(Feature_ID=feature_data$Feature_ID, RT_err=abs(reference$rt[k]-feature_data$rt), stringsAsFactors = F)
    res[[k]]$RT_flag <- res[[k]]$RT_err < 1
    res[[k]] <- res[[k]][res[[k]]$RT_flag, ]
  }

  return(res)
}
