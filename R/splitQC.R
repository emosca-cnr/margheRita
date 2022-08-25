#' Split dataframe and metadata in two: one for samples and one for QC

#' @param mRList mRList object
#' @export
#' @return


splitQC<-function(mRList){

  idx_QC <- mRList$sample_ann$class=="QC"

  mRList$QC <- mRList$data[, idx_QC]
  mRList$QC_ann <- mRList$sample_ann[idx_QC, ]

  mRList$data <- mRList$data[, !idx_QC]
  mRList$sample_ann <- mRList$sample_ann[!idx_QC, ]


  return(mRList)

}

