#' Merge feature files with statistics
#' @description Merge feature files with results of the statistics, and export the results
#' @param mRList mRList object
#' @param test_method name of an element of mRList object containing results of statistical test or a custom data frame with Feature_ID as row names (e.g.: "anova").
#' @param test_value value to use for significance. e.g.: "p" or "q"
#' @param cutoff_value value below which a feature is considered significant.
#' @param dirout output directory
#' @return the data.frame metab_stat, metab_stat_w_features, metab_stat_w_features_known, metab_stat_w_features_only_sign, metab_stat_w_features_known_only_sign will be both saved in the global environemtn and exported as txt tab delimited files.
#' @export


merge_stats_w_features <- function(mRList = NULL, test_method = "anova", test_value = "q", cutoff_value = 0.05, dirout=NULL){
  
  stopifnot(!is.null(test_method))
  
  if (!is.null(dirout)) {
    dir.create(dirout, showWarnings = F)
  }else{
    dirout <- getwd()
  }
  
  if(!is.data.frame(test_method)){
    tab_test_method <- mRList[[test_method]]
  }
  
  #metab_ann contains the merged features...
  #merge annotation with statistical analysis
  metab_stat <- merge(mRList$metab_ann, tab_test_method, by.x = "Feature_ID", by.y=0)
  metab_stat <- metab_stat[, ! colnames(metab_stat) %in% c("rt", "mz", "MS1.isotopic.spectrum", "MS_MS_spectrum", "quality", "reference")]
  
  ans <- list(metab_stat=metab_stat)
  rm(metab_stat)
  
  ## merge metab_stat con i dati delle features data
  ans$metab_stat_feat <- merge(ans$metab_stat, mRList[["data"]], by.x="Feature_ID", by.y=0, sort=F)

  ## Only known
  ans$metab_stat_feat_known <- ans$metab_stat_feat[ans$metab_stat_feat$Feature_ID %in% mRList$data_ann$Feature_ID, ]
  
  
  ### SIGNIFICANT FEATURES
  sign_feat_featID <- select_sign_features(mRList = mRList, test_method = test_method, test_value = test_value, cutoff_value = cutoff_value, feature_id = "Feature_ID")
  
  
  ## Only significant
  ans$metab_stat_feat_sig <- ans$metab_stat_feat[ans$metab_stat_feat$Feature_ID %in% sign_feat_featID, ]
  
  ## Only known & significant
  ans$metab_stat_feat_known_sig <- ans$metab_stat_feat_known[ans$metab_stat_feat_known$Feature_ID %in% ans$metab_stat_feat_known$Feature_ID, ]
  
  ### write to files
  wb <- createWorkbook()
  for(i in 1:length(ans)){
    #write.table(ans[[i]], file = paste0(names(ans)[i], ".2.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
    addWorksheet(wb, names(ans)[i])
    writeDataTable(wb, names(ans)[i], ans[[i]])
  }
  saveWorkbook(wb, file = file.path(dirout, "features_statistics_annotation.xlsx"), overwrite = T)
  
  return(ans)
  
}


