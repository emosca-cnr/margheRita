#' calculate reference
#' @export
#'
calc_reference <- function(m_list, sample_col="type", sample_type="QC", approach=c("median", "mean")){

  approach <- approach[1]

  idx_samples <- m_list$sample_ann[, colnames(m_list$sample_ann)==sample_col]
  idx_samples <- idx_samples==sample_type

  if(approach == "median"){
   m_list$metab_ann$reference <- apply(m_list$data[, idx_samples], 1, median)
  }
  if(approach == "mean"){
    m_list$metab_ann$reference <- RowMeans(m_list$data[, idx_samples])
  }

  return(m_list)
}
