#' Calculate of coefficient of variation for each metabolite
#'
#' @description Calculate the coefficient of variance (CV) for samples and QC (quality control) to assess the reliability of the experiment, as:
#' CV = (SD - mean) * 100.
#' CV samples and CV QC are compared, calculating the ratio between CV samples and CV QC.
#' The threshold of the ratio has to be decided by the user in ratioCV.
#' Only metabolites with ratio >= ratioCV are used for the analysis.
#'
#' @param mRList mRList object
#' @param ratioCV User set the value of the ratio between CV samples/ CC QC
#' @export
#' @return filtered mRList object
#' @param dirout output directory
#' @examples
#' ##library(dataset.margheRita)
#' ##dataset(norm_pos)
#' mRList<-CV(mRList, dirout, ratioCV=2)
#'

CV_ratio <- function(mRList, dirout="./", ratioCV) {

  dir.create(dirout)

  mean_QC <- apply(mRList$QC, MARGIN=1, mean)
  sd_QC <- apply(mRList$QC, MARGIN = 1, sd)
  CV_QC <- (sd_QC / mean_QC) * 100
  utils::write.csv(CV_QC, file = paste0(dirout,"/CV_QC.csv"))


  mean_Samples <- apply(mRList$data, MARGIN=1, mean)
  sd_Samples <- apply(mRList$data, MARGIN = 1, sd)
  CV_Samples <- (sd_Samples / mean_Samples) * 100
  utils::write.csv(CV_Samples, file = paste0(dirout,"/CV_Sample.csv"))

  CV_all <- cbind(CV_QC, CV_Samples)
  colnames(CV_all) = c("CV_QC", "CV_Samples")
  utils::write.csv(CV_all, paste0(dirout, "/CV_all.csv"))

  cat("Summary of CV_samples / CV_QC:\n")
  print(summary(CV_all[, 1]/ CV_all[, 2]))

  ratio<-CV_Samples/CV_QC
  ratio<-as.data.frame(ratio)
  rownames(ratio)<-mRList$metab_ann$Feature_ID

  CV_all <- cbind(CV_Samples,CV_QC,ratio)
  colnames(CV_all) = c( "CV_Samples", "CV_QC", "ratio")
  utils::write.csv(CV_all, paste0(dirout, "/CV_all.csv"))

  idx_keep <- CV_Samples/CV_QC >= ratioCV
  cat("# Metabolites with appropriate CV", sum(idx_keep), "/", nrow(CV_all), "\n")

  mRList$data <- mRList$data[idx_keep, ] #mRList$data cleaned
  mRList$metab_ann <- mRList$metab_ann[idx_keep, ] #mRList$data cleaned

  mRList$QC <- mRList$QC[idx_keep, ] #mRList$data cleaned

  return(mRList)
}


