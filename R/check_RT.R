#' Check the retention time
#'
#' @param reference data frame with retention times (rt) of reference metabolites.
#' @param feature_data data frame with observed rt.
#' @param rt_err_thr threshold for acceptable rt.
#' @param rt_best_thr threshold for super rt.
#' @param filter whether to consider only features with RT error <= rt_err_thr
#'
#' @return A list of where each element is named with a reference library ID and contains a data frame of features
#'


check_RT <- function(reference=NULL, feature_data=NULL, rt_err_thr = 1, rt_best_thr=0.5, filter=FALSE){
  
  RT <- vector("list", nrow(reference))
  names(RT) <- reference$ID
  
  for (k in 1:dim(reference)[1]){
    
    RT[[k]] <- data.frame(Feature_ID= feature_data$Feature_ID, rt=feature_data$rt, RT_err=abs(reference$rt[k]- feature_data$rt), stringsAsFactors = F)
    RT[[k]]$RT_flag <- RT[[k]]$RT_err <= rt_err_thr
    
    RT[[k]]$RT_class <- "unacceptable"
    RT[[k]]$RT_class[ RT[[k]]$RT_flag] <- "acceptable"
    RT[[k]]$RT_class[ RT[[k]]$RT_err < rt_best_thr] <- "super"
    
    if(filter){
      RT[[k]] <- RT[[k]][RT[[k]]$RT_flag, ]
    }

  }
  
  return(RT)
}
