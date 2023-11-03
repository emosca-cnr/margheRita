#' Check the retention time
#'
#' @description `check_RT` it is making a list of all library metabolites and assigning proper sample ID to each based on desired retention time error.
#'
#' @param reference A list of library contain retention time with specific ID.
#' @param feature_data A list of sample data contain retention time with specific ID.
#' @param rt_err_thr A number that specify the desired absolute difference of retention time between sample and library metabolite.
#' @param rt_best_thr Features with RT error < rt_best_thr  will be labelled as "super", those with rt_best_thr < RT error <= rt_err_thr as "acceptable", the remaining as "unacceptable".
#' @param filter whether to consider only features with RT error <= rt_err_thr
#'
#' @return A list of library ID each contain a data frame of sample ID with retention time in range of rt_err_thr
#' @export
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
