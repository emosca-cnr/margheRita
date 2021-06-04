#' calculate reference
#' @export
#'
calc_reference <- function(m_list, sample_col="class", sample_class="QC", approach=c("median", "mean")){

  approach <- approach[1]

  if(sample_class=="QC" & !is.null(m_list$QC)){
    cat("Using QC...\n")
    X <- m_list$QC

  }else{

    cat("Using ", sample_col, "of class ", sample_class, "\n")
    idx_samples <- m_list$sample_ann[, colnames(m_list$sample_ann)==sample_col]
    idx_samples <- idx_samples==sample_class

    X <- m_list$data[, idx_samples]

  }
  if(approach == "median"){
    m_list$metab_ann$reference <- apply(X, 1, median)
  }
  if(approach == "mean"){
    m_list$metab_ann$reference <- RowMeans(X)
  }

  return(m_list)
}
