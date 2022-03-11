#' calculate reference
#' @export
#'
calc_reference <- function(mRList, sample_col="class", sample_class="QC", approach=c("median", "mean")){

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
    mRList$metab_ann$reference <- RowMeans(X)
  }

  return(mRList)
}
