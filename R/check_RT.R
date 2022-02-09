#' check_RT
#'
#'
#' @export


check_RT <- function(feature_data, reference, rt_err_thr = rt_err_thr){

  RT <- vector("list", nrow(reference))
  names(RT) <- reference$Name

  for (k in 1:dim(reference)[1]){

    #idx_ok <- which((feature_data$rt > (reference$rt[k] - rt_window)) & (feature_data$rt < (reference$rt[k] + rt_window)))

    RT[[k]] <- data.frame(Feature_ID= feature_data$Feature_ID, RT_err=abs(reference$rt[k]- feature_data$rt), stringsAsFactors = F)
    RT[[k]]$RT_flag <- RT[[k]]$RT_err < rt_err_thr
    RT[[k]] <- RT[[k]][RT[[k]]$RT_flag, ]
  }

  #filter the RT by deleting the empty data.frame;
  #RT = RT[sapply(RT, function(x) dim(x)[1]) > 0]


  return(RT)
}
