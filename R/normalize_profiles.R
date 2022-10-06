#' Normalize profiles
#' @param mRList mRList object
#' @param method "log": log2(x+1); "reference": divide each sample by its reference; "pqn": probailistic quotient normalization;
#' @export
#' @return mRList object with normalized data

normalize_profiles <- function(mRList, method=c("log", "reference", "pqn")){

  if(method == "log"){
    mRList$data <- log2(mRList$data+1)
  }

  if(method == "reference"){
    mRList$data <- t(t(mRList$data) / mRList$metab_ann$reference)
  }

  if(method == "pqn"){
    mRList$data <- apply(mRList$data, 2, function(x) pqn(x, xref = mRList$metab_ann$reference)$y)
  }

  return(mRList)

}
