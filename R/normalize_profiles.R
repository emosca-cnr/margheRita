#' Normalize profiles
#' @param mRList mRList object
#' @param method "log": log2(x+1); "reference": divide each sample by its reference; "pqn": probailistic quotient normalization;
#' @export
#' @return mRList object with normalized data

normalize_profiles <- function(mRList=NULL, method=c("log", "reference", "pqn")){

  if(method == "log"){
    cat("log normalization\n")
    mRList$data <- log2(mRList$data+1)
  }

  if(method == "reference"){
    if(length(mRList$metab_ann$reference)==0){
      cat("No reference profile found, using calc_reference() function...\n")
      mRList <- calc_reference(mRList)
    }
    cat("normalizing to a reference\n")
    mRList$data <- t(t(mRList$data) / mRList$metab_ann$reference)
  }

  if(method == "pqn"){
    cat("PQN normalization\n")
    if(length(mRList$metab_ann$reference)==0){
      cat("No reference profile found, using calc_reference() function...\n")
      mRList <- calc_reference(mRList)
    }
    mRList$data <- apply(mRList$data, 2, function(x) pqn(x, xref = mRList$metab_ann$reference)$y)
  }

  return(mRList)

}
