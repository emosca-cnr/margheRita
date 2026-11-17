#' Coefficient of variation filter
#' Calculate the ratio between CVs (coefficient of variation) of samples and QCs to assess the reliability of the experiment. Only features with ratio >= ratioCV are kept. It produces a text file with in dirout.
#'
#' @param mRList mRList object
#' @param ratioCV only metabolites with CV ratio (sample/QC) higher than ratioCV will be kept
#' @param dirout output directory
#' @importFrom utils write.csv
#' @importFrom stats sd
#' @return mRList object
#' @export

CV_ratio <- function(mRList=NULL, dirout=NULL, ratioCV=1) {

  if (!is.null(dirout)) {
    dir.create(dirout, showWarnings = F)
  }else{
    dirout <- getwd()
  }


  if(!identical(rownames(mRList$data), rownames(mRList$QC))){
    stop("Rownames of data and QC are not identical.\n")
  }

  mean_QC <- apply(mRList$QC, MARGIN=1, mean)
  sd_QC <- apply(mRList$QC, MARGIN = 1, sd)
  CV_QC <- (sd_QC / mean_QC)
  
  #if some feature is equal in all WC samnples  CV_QC is null: it is set to the min value of the dataset
  if(any(CV_QC==0)){
    CV_QC[CV_QC==0] <- min(CV_QC[CV_QC!=0])
  }
  
  mean_Samples <- apply(mRList$data, MARGIN=1, mean)
  sd_Samples <- apply(mRList$data, MARGIN = 1, sd)
  CV_Samples <- (sd_Samples / mean_Samples)

  ratio <- CV_Samples / CV_QC

  CV_all <- data.frame(Feature_ID=names(CV_Samples), Samples=CV_Samples, QC=CV_QC, CVr=ratio, stringsAsFactors = F)
  write.csv(CV_all, file.path(dirout, "CV_all.csv"), row.names = F)

  cat("Summary of CV ratio (samples / QC):\n")
  print(summary(ratio))

  idx_keep <- ratio > ratioCV
  cat("# Features with appropriate CV", sum(idx_keep), "/", nrow(CV_all), "\n")

  mRList$data <- mRList$data[idx_keep, ] #mRList$data cleaned
  mRList$metab_ann <- mRList$metab_ann[idx_keep, ] #mRList$data cleaned

  mRList$QC <- mRList$QC[idx_keep, ] #mRList$data cleaned

  return(mRList)
}


