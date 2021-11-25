#' check_RT
#'
#'
#' @export


check_RT <- function(feature_data, reference){

  RT <- vector("list", nrow(reference))
  names(RT) <- reference$Name

  for (k in 1:dim(reference)[1]){

    idx_ok <- which((feature_data$rt > (reference$rt[k] - 2)) & (feature_data$rt < (reference$rt[k] +2)))

    RT[[k]] <- data.frame(Feature_ID=feature_data$Feature_ID, RT_err=abs(reference$rt[k]-feature_data$rt), stringsAsFactors = F)
    RT[[k]]$RT_flag <- RT[[k]]$RT_err < 2
    RT[[k]] <- RT[[k]][RT[[k]]$RT_flag, ]
  }

  return(RT)
}
