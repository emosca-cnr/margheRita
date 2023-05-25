#' Calculate a reference profile
#'
#' The reference profile is calculated using the samples annotated as "sample_class" in the sample annotation column "sample_col".
#' Launch function before normalize_profiles function
#' @param mRList mRList object
#' @param sample_col the column of mRList$sample_ann where the "sample_class" value is found
#' @param sample_class label to identify the samples to be used for the reference profile calculation
#' @param approach The type of average: "median" or "mean".
#'
#' @export
#' @return mRList object with mRList$reference
#' @examples
#' ##library(dataset.margheRita)
#' ##dataset(norm_pos)
#' ## mRList<-calc_reference(mRList=mRList, sample_col="class", sample_class="QC", approach="mean")


calc_reference <- function(mRList=NULL, sample_col="class", sample_class="QC", approach=c("median", "mean")){

  approach <- approach[1]

  if(sample_class=="QC" & !is.null(mRList$QC)){
    cat("Using QC...\n")
    X <- mRList$QC

  }else{

    cat("Using ", sample_col, "of class ", sample_class, "\n")
    idx_samples <- mRList$sample_ann[, colnames(mRList$sample_ann)==sample_col]
    idx_samples <- idx_samples==sample_class

    X <- mRList$data[, idx_samples]

  }
  if(approach == "median"){
    mRList$metab_ann$reference <- apply(X, 1, median)
  }
  if(approach == "mean"){
    mRList$metab_ann$reference <- rowMeans(X)
  }

  return(mRList)
}
