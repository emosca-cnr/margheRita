#' Check the retention time
#'
#' @description `check_RT` it is making a list of all library metabolites and assigning proper sample ID to each based on desired retention time error.
#'
#' @param reference A list of library contain retention time with specific ID.
#' @param feature_data A list of sample data contain retention time with specific ID.
#' @param rt_err_thr A number that specify the desired absolute difference of retention time between sample and library metabolite.
#'
#' @return A list of library ID each contain a data frame of sample ID with retention time in range of rt_err_thr
#' @export
#'


check_RT <- function(reference=NULL, feature_data=NULL, rt_err_thr = 1){

  RT <- vector("list", nrow(reference))
  names(RT) <- reference$ID

  for (k in 1:dim(reference)[1]){

    RT[[k]] <- data.frame(Feature_ID= feature_data$Feature_ID, RT_err=abs(reference$rt[k]- feature_data$rt), stringsAsFactors = F)
    RT[[k]]$RT_flag <- RT[[k]]$RT_err < rt_err_thr
    RT[[k]] <- RT[[k]][RT[[k]]$RT_flag, ]
  }

  return(RT)
}
