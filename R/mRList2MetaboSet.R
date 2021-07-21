#' mRList2MetaboSet
#'
#' @import notame
#' @export

mRList2MetaboSet <- function(mRList){

  ### 1.1 #### create a notame MetaboSet 0bject
  feature_data <- m_list_init$metab_ann
  colnames(feature_data)[1] <- "Feature_ID"
  colnames(feature_data)[2] <- "rt"
  colnames(feature_data)[3] <- "Average mz"
  feature_data <- cbind(feature_data, Split="HILIC_pos")
  colnames(feature_data)[1] <- "Feature_ID"

  pheno_data <- rbind(m_list_init$sample_ann, m_list_init$QC_ann)
  colnames(pheno_data)[1] <- "Sample_ID"
  colnames(pheno_data)[3] <- "Injection_order"
  colnames(pheno_data)[5] <- "QC"

  mset <- notame::construct_metabosets(exprs = cbind(m_list_init$data, m_list_init$QC), pheno_data = pheno_data, feature_data = feature_data, group_col = "class")

  return(mset)

}
