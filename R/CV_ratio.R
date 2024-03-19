#' Coefficient of variation filter
#'
#' @description Calculate the coefficient of variance (CV) for samples and QC (quality control) to assess the reliability of the experiment, as:
#' CV = (SD - mean) * 100.
#' CV samples and CV QC are compared, calculating the ratio between CV samples and CV QC.
#' The threshold of the ratio has to be decided by the user in ratioCV.
#' Only metabolites with ratio >= ratioCV are used for the analysis.
#'
#' @param mRList mRList object
#' @param ratioCV only metabolites with CV ratio (sample/QC) higher than ratioCV will be kept
#' @export
#' @return filtered mRList object
#' @param dirout output directory
#' @importFrom utils write.csv

CV_ratio <- function(mRList=NULL, dirout="./", ratioCV=1) {

  if (dirout != "./") {
    dir.create(dirout)
  }


  if(!identical(rownames(mRList$data), rownames(mRList$QC))){
    stop("Rownames of data and QC are not identical.\n")
  }

  mean_QC <- apply(mRList$QC, MARGIN=1, mean)
  sd_QC <- apply(mRList$QC, MARGIN = 1, sd)
  CV_QC <- (sd_QC / mean_QC)

  mean_Samples <- apply(mRList$data, MARGIN=1, mean)
  sd_Samples <- apply(mRList$data, MARGIN = 1, sd)
  CV_Samples <- (sd_Samples / mean_Samples)

  ratio <- CV_Samples / CV_QC

  CV_all <- data.frame(Feature_ID=names(CV_Samples), Samples=CV_Samples, QC=CV_QC, CVr=ratio, stringsAsFactors = F)
  write.csv(CV_all, file.path(dirout, "CV_all.csv"), row.names = F)

  cat("Summary of CV ratio (samples / QC):\n")
  print(summary(ratio))

  idx_keep <- ratio > ratioCV
  cat("# Metabolites with appropriate CV", sum(idx_keep), "/", nrow(CV_all), "\n")

  mRList$data <- mRList$data[idx_keep, ] #mRList$data cleaned
  mRList$metab_ann <- mRList$metab_ann[idx_keep, ] #mRList$data cleaned

  mRList$QC <- mRList$QC[idx_keep, ] #mRList$data cleaned

  return(mRList)
}


