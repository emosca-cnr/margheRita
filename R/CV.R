#' Calculation of coefficient of variance for each metabolite
#' @description Function Calculates coefficient of variance (CV) for samples and QC (quality control) to assess the reliability of the experiment.
#' @description CV as =(stdev-mean)*100, is calculated for each metabolite for QC and samples. If CV QC>CV sample metabolite is deleted.
#' @description The resulting dataframe and metadata contain only reliable metabolites.
#' @param mRList mRList has to be split into df_QC (dataframe of quality control) and df_samples (dataframe of samples)
#' @export
#' @return mRList cleaned with only reliable metabolites



#apply CV QC vs CV of samples: eliminate when CV of QC >CV of samples for each metabolites

CV <- function(mRList, dirout) {

  dir.create(dirout)

  mean_QC <- apply(mRList$QC, MARGIN=1, mean)
  sd_QC <- apply(mRList$QC, MARGIN = 1, sd)
  CV_QC <- (sd_QC / mean_QC) * 100
  utils::write.csv(CV_QC, file = paste0(dirout,"/CV_QC.csv"))


  mean_Samples <- apply(mRList$data, MARGIN=1, mean)
  sd_Samples <- apply(mRList$data,MARGIN = 1,sd)
  CV_Samples <- (sd_Samples / mean_Samples) * 100
  utils::write.csv(CV_Samples, file = paste0(dirout,"/CV_Sample.csv"))

  CV_all <- cbind(CV_QC, CV_Samples)
  colnames(CV_all) = c("CV_QC", "CV_Samples")
  utils::write.csv(CV_all, paste0(dirout, "/CV_all.csv"))

  idx_keep <- CV_QC < CV_Samples
  cat("# Metabolites with appropriate CV", sum(idx_keep), "\n")

  mRList$data <- mRList$data[idx_keep, ] #mRList$data cleaned
  mRList$metab_ann <- mRList$metab_ann[idx_keep, ] #mRList$data cleaned

  mRList$QC <- mRList$QC[idx_keep, ] #mRList$data cleaned

  return(mRList)
}


