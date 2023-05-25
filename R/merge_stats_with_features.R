#' Get significant features
#' Merge mRList$metab_ann with a data.frame with statistical analysis
#' @param mRList mRList object
#' @param feature_stats name of an element of mRList object containing results of statistical test or a custom data frame with Feature_ID as row names.
#' @export

merge_stats_with_features <- function(mRList = NULL, feature_stats=NULL){
  
  if(is.null(feature_stats)){
    stop("feature_stats argument can not be NULL.\n")
  }
  if(!is.data.frame(feature_stats)){
    feature_stats <- mRList[[feature_stats]]
  }
  
  ans <- merge(mRList_norm_bio$metab_ann, feature_stats, by.x = "Feature_ID", by.y=0)
  ans <- ans[, ! colnames(ans) %in% c("rt", "mz", "MS1.isotopic.spectrum", "MS_MS_spectrum", "quality", "reference")]
  
  return(ans)
  
}