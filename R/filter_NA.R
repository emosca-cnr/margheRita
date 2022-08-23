#' filter NA
#' @param mRList mRList object
#' @param min_metab_in_sample min number of metabolites in a sample
#' @param min_sample_with_metab min number of samples in which a metabolite must appear
#' @export
#' @return filtered mRList object


filter_NA <- function(mRList, min_metab_in_sample=100, min_sample_with_metab=10){

  idx_keep <- colSums(!is.na(mRList$data)) >= min_metab_in_sample
  cat("# Samples with enough metabolites", sum(idx_keep), "/", ncol(mRList$data), "\n")
  mRList$data <- mRList$data[, idx_keep]
  mRList$sample_ann <- mRList$sample_ann[idx_keep, ]

  idx_keep <- rowSums(!is.na(mRList$data)) >= min_sample_with_metab
  cat("# Metabolites occurring in enough samples", sum(idx_keep), "/", nrow(mRList$data), "\n")

  mRList$data <- mRList$data[idx_keep, ]
  mRList$metab_ann <- mRList$metab_ann[idx_keep, ]

  return(mRList)
}
