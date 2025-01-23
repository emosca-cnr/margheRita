#' Merge feature files with statistics
#' @description Merge feature files with results of the statistics, and export the results
#' @param mRList mRList object
#' @param test_method name of an element of mRList object containing results of statistical test or a custom data frame with Feature_ID as row names (e.g.: "anova").
#' @param test_value value to use for significance. e.g.: "p" or "q"
#' @param cutoff_value value below which a feature is considered significant.
#' @return the data.frame metab_stat, metab_stat_w_features, metab_stat_w_features_known, metab_stat_w_features_only_sign, metab_stat_w_features_known_only_sign will be both saved in the global environemtn and exported as txt tab delimited files.
#' @importFrom dplyr left_join
#' @export


merge_stats_w_features <- function(mRList = NULL, test_method = "anova", test_value = "q", cutoff_value = 0.05){
  
  if(is.null(test_method)){
    stop("feature_stats argument can not be NULL.\n")
  }
  if(!is.data.frame(test_method)){
    tab_test_method <- mRList[[test_method]]
  }
  
  ans <- merge(mRList$metab_ann, tab_test_method, by.x = "Feature_ID", by.y=0)
  ans <- ans[, ! colnames(ans) %in% c("rt", "mz", "MS1.isotopic.spectrum", "MS_MS_spectrum", "quality", "reference")]
  
  data_annotation <- as_tibble(mRList[["metabolite_identification"]][["associations"]])
  # order data annotation not only per Level, but also for the other characteristic, to better prioritize:
  data_annotation$Level <- as.numeric(data_annotation$Level)
  data_annotation$mass_status <- factor(data_annotation$mass_status, levels = c("suffer", "acceptable", "super"))
  data_annotation$peaks_found_ppm_RI <- as.numeric(data_annotation$peaks_found_ppm_RI)
  data_annotation$matched_peaks_ratio <- as.numeric(data_annotation$matched_peaks_ratio)
  data_annotation$RT_class <- factor(data_annotation$RT_class, levels = c("unacceptable", "acceptable", "super"))
  data_annotation$ppm_error <- as.numeric(data_annotation$ppm_error)
  data_annotation$RT_err <- as.numeric(data_annotation$RT_err)
  
  data_annotation_ordered <- arrange(data_annotation, Level, desc(mass_status), desc(peaks_found_ppm_RI), desc(matched_peaks_ratio), desc(RT_class), ppm_error, RT_err)
  
  
  
  metab_stat <- left_join(x = ans, y = data_annotation_ordered, by = "Feature_ID", suffix = c("", "_bis"), multiple = "first")
  write.table(metab_stat,
              file = "metab_stat.txt",
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE)
  assign("metab_stat", metab_stat, envir = .GlobalEnv)
  
  
  ## merge metab_stat con i dati delle features
  
  
  data_intensities <- cbind(data.frame(Feature_ID = row.names(mRList[["data"]])),
                            as.data.frame(mRList[["data"]]))
  
  
  metab_stat_w_features <- left_join(x = metab_stat, y = data_intensities, by = "Feature_ID")
  write.table(metab_stat_w_features,
              file = "metab_stat_w_features.txt",
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE)
  assign("metab_stat_w_features", metab_stat_w_features, envir = .GlobalEnv)
  
  
  
  
  ## uguale ma solo quelle known:
  
  
  data_intensities_known <- cbind(data.frame(Name = row.names(mRList[["data_ann"]])),
                                  as.data.frame(mRList[["data_ann"]]))
  colnames(data_intensities_known)[duplicated(colnames(data_intensities_known)) & colnames(data_intensities_known)== "Name"] <- "Name_bis"
  
  
  metab_stat$Level <- as.numeric(metab_stat$Level)
  metab_stat[,test_value] <- as.numeric(pull(metab_stat, test_value))
  
  metab_stat <- arrange(metab_stat, Level, !!sym(test_value))
  metab_stat_w_features_known <- left_join(x = data_intensities_known, y = metab_stat, by = "Feature_ID", multiple = "first")
  write.table(metab_stat_w_features_known,
              file = "metab_stat_w_features_known.txt",
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE)
  assign("metab_stat_w_features_known", metab_stat_w_features_known, envir = .GlobalEnv)
  
  
  ## uguale ma solo quelle significative:
  
  sign_feat_featID <- select_sign_features(mRList = mRList, test_method = test_method, test_value = test_value, cutoff_value = cutoff_value, feature_id = "Feature_ID")
  
  data_intensities_only_sign <- data_intensities[which(data_intensities$Feature_ID %in% sign_feat_featID),]
  
  
  metab_stat_w_features_only_sign <- left_join(x = data_intensities_only_sign , y = metab_stat, by = "Feature_ID", multiple = "first")
  write.table(metab_stat_w_features_only_sign,
              file = "metab_stat_w_features_only_sign.txt",
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE)
  assign("metab_stat_w_features_only_sign", metab_stat_w_features_only_sign, envir = .GlobalEnv)
  
  
  
  ## uguale ma solo quelle known e significative:
  
  
  data_intensities_known_only_sign <- data_intensities_known[which(data_intensities_known$Feature_ID %in% data_intensities_only_sign$Feature_ID),]
  
  
  metab_stat_w_features_known_only_sign <- left_join(x = data_intensities_known_only_sign , y = metab_stat, by = "Feature_ID", multiple = "first")
  write.table(metab_stat_w_features_known_only_sign,
              file = "metab_stat_w_features_known_only_sign.txt",
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE)
  assign("metab_stat_w_features_known_only_sign", metab_stat_w_features_known_only_sign, envir = .GlobalEnv)
  
  
}


