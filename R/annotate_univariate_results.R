#' Get significant features
#' @description Merge mRList$metab_ann with feature_stats
#' @param mRList mRList object
#' @param feature_stats name of an element of mRList object containing results of statistical test or a custom data frame with Feature_ID as row names.
#' @param add.metabolite.levels if TRUE, adds mRList$data_ann, if this exists.
#' @param dirout output directory
#' @return a data.frame
#' @export

annotate_univariate_results <- function(mRList = NULL, feature_stats=NULL, add.metabolite.levels=TRUE, dirout=NULL){
  
  if(is.null(feature_stats)){
    stop("feature_stats argument can not be NULL.\n")
  }
  
  if (!is.null(dirout)) {
    dir.create(dirout, showWarnings = F)
  }else{
    dirout <- getwd()
  }
  
  
  if(!is.data.frame(feature_stats)){
    stopifnot(length(mRList[[feature_stats]])>0)
    feature_stats <- cbind(Feature_ID=rownames(mRList[[feature_stats]]), mRList[[feature_stats]])
  }
  
  if(length(mRList$metabolite_identification) > 0){
    
    if(add.metabolite.levels){
      ans <- merge(mRList$data_ann, feature_stats, by = "Feature_ID")
      ans <- merge(ans, mRList$metabolite_identification$associations, by = c("Feature_ID", "Name"))
    }else{
      ans <- merge(feature_stats, mRList$metabolite_identification$associations, by = "Feature_ID")
    }
    
  }else{
    
    ans <- feature_stats
    if(add.metabolite.levels){
      ans <- merge(cbind(Feature_ID=rownames(mRList$data), mRList$data), ans, by = "Feature_ID")
    }
    ans <- merge(ans, mRList$metab_ann, by = "Feature_ID")
    ans <- ans[, ! colnames(ans) %in% c("rt", "mz", "MS1.isotopic.spectrum", "MS_MS_spectrum", "quality", "reference")]
    
  }
  
  write.csv(ans, file = file.path(dirout, "data_stats_ann.csv"), row.names = F)
  
  return(ans)
  
}