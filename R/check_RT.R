#' check_RT
#'
#'
#' @export


check_RT <- function(feature_data, reference){

  feature_data$PEPMASS <- as.numeric(feature_data$PEPMASS)
  feature_data$RTINMINUTES <- as.numeric(feature_data$RTINMINUTES)

  res <- vector("list", nrow(reference))
  names(res) <- reference$`Metabolite name`

  for (k in 1:dim(reference)[1]){

    idx_ok <- which((feature_data$RTINMINUTES > (reference$`Average Rt(min)`[k] - 1)) & (feature_data$RTINMINUTES < (reference$`Average Rt(min)`[k] +1)))

    res[[k]] <- data.frame(ID=feature_data$ID, RT_err=abs(reference$`Average Rt(min)`[k]-feature_data$RTINMINUTES), stringsAsFactors = F)
    res[[k]]$RT_flag <- res[[k]]$RT_err < 1
    res[[k]] <- res[[k]][res[[k]]$RT_flag, ]
  }

  return(res)
}
