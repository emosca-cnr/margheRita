#' filter NA
#' @param m_list margheRita list
#' @param min_metab_in_sample min number of metabolites in a sample
#' @param min_sample_with_metab min number of samples in which a metabolite must appear
#' @export
filter_NA <- function(m_list, min_metab_in_sample=100, min_sample_with_metab=10){

  idx_keep <- colSums(!is.na(m_list$df)) >= min_metab_in_sample
  m_list$df <- m_list$df[, idx_keep]
  m_list$sample_ann <- m_list$sample_ann[, idx_keep]

  idx_keep <- rowSums(!is.na(m_list$df)) >= min_sample_with_metab
  m_list$df <- m_list$df[idx_keep, ]
  m_list$metab_ann <- m_list$metab_ann[idx_keep, ]

  return(m_list)
}
