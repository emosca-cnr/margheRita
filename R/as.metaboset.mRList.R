#' Transforms mRList into mSet object
#' @param mRList mRList object
#' @importFrom notame construct_metabosets
#' @export
#' @return mset object

as.metaboset.mRList <- function(mRList){

  ### 1.1 #### create a notame MetaboSet 0bject
  feature_data <- mRList$metab_ann
  colnames(feature_data)[3] <- "Average mz"
  feature_data <- cbind(feature_data, Split="HILIC_pos")

  pheno_data <- rbind(mRList$sample_ann, mRList$QC_ann)
  colnames(pheno_data)[1] <- "Sample_ID"
  colnames(pheno_data)[colnames(pheno_data) == "injection_order"] <- "Injection_order"
  pheno_data <- cbind(pheno_data, QC=pheno_data$class)

  mset <- notame::construct_metabosets(exprs = as.matrix(cbind(mRList$data, mRList$QC)), pheno_data = pheno_data, feature_data = feature_data, group_col = "class")

  return(mset)

}
