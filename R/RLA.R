#' Relative Log Abudance
#' @param mRList mRList object
#' @param include_QC whether to include or not the QC samples
#' @param logged are the input data on log-scale or not?
#' @param robust whether to use the median or not
#' @importFrom stats median
#' @return mRList object with element "RLA"
#' @export
#' 
RLA <- function(mRList=NULL, include_QC=FALSE, logged=TRUE, robust=TRUE){


  if(include_QC){
    X_data <- cbind(mRList$data, mRList$QC)
    X_sample_ann <- rbind(mRList$sample_ann, mRList$QC_ann)
  }else{
    X_data <- mRList$data
    X_sample_ann <- mRList$sample_ann
  }

  if(!logged){
    X_data <- log2(X_data + 1)
  }

  if(robust){
    X_data <- t(apply(X_data, 1, function(y) y - median(y)))
  }else{
    X_data <- t(apply(X_data, 1, function(y) y - mean(y)))
  }

  mRList$RLA <- X_data

  return(mRList)

}
