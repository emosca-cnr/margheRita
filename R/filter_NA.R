#' filter NA
#' @param mRList margheRita list
#' @param min_metab_in_sample min number of metabolites in a sample
#' @param min_sample_with_metab min number of samples in which a metabolite must appear
#' @export
filter_NA <- function(mRList, min_metab_in_sample=100, min_sample_with_metab=10){

  idx_keep <- colSums(!is.na(mRList$data)) >= min_metab_in_sample
  cat("Samples with enough metabolites\n")
  print(table(idx_keep))
  mRList$data <- mRList$data[, idx_keep]
  mRList$sample_ann <- mRList$sample_ann[idx_keep, ]

  idx_keep <- rowSums(!is.na(mRList$data)) >= min_sample_with_metab
  cat("Metabolites occurring in enough samples\n")
  print(table(idx_keep))
  mRList$data <- mRList$data[idx_keep, ]
  mRList$metab_ann <- mRList$metab_ann[idx_keep, ]

  return(mRList)
}
