#' filter NA
#' @param mRList mRList object
#' @param min_metab_in_sample min number of metabolites in a sample
#' @param min_sample_with_metab min number of samples in which a metabolite must appear
#' @param na_value value that indicate missing values
#' @export
#' @return filtered mRList object


filter_NA <- function(mRList=NULL, min_metab_in_sample=100, min_sample_with_metab=3, na_value="NA"){

  if(na_value != "NA"){
    cat("setting", na_value, "to NA\n")
    mRList$data[mRList$data == na_value] <- NA
  }

  idx_keep <- colSums(!is.na(mRList$data)) >= min_metab_in_sample
  cat("# Samples with >=", min_metab_in_sample, "metabolites", sum(idx_keep), "/", ncol(mRList$data), "\n")
  mRList$data <- mRList$data[, idx_keep]
  mRList$sample_ann <- mRList$sample_ann[idx_keep, ]

  idx_keep <- rowSums(!is.na(mRList$data)) >= min_sample_with_metab
  cat("# Features occurring in >= ", min_sample_with_metab, "samples", sum(idx_keep), "/", nrow(mRList$data), "\n")

  mRList$data <- mRList$data[idx_keep, ]
  mRList$metab_ann <- mRList$metab_ann[idx_keep, ]

  return(mRList)
}
