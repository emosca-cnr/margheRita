#' Transforms mRList into MetaboSet object
#' @param mRList mRList object
#' @export
#' @return metabosets
#' @importFrom notame construct_metabosets
#' @importFrom Biobase AnnotatedDataFrame fData

as.metaboset <- function(mRList=NULL){
  
  ### 1.1 #### create a notame MetaboSet 0bject
  feature_data <- mRList$metab_ann
  colnames(feature_data)[colnames(feature_data) == "mz"] <- "Average mz"
  feature_data$Split <- FALSE
  
  pheno_data <- mRList$sample_ann
  exprs_data <- mRList$data
  
  if(length(mRList$QC_ann) > 0){
    pheno_data <- rbind(mRList$sample_ann, mRList$QC_ann)
    pheno_data <- cbind(pheno_data, QC=pheno_data$class)
    exprs_data <- cbind(exprs_data, mRList$QC)
  }
  colnames(pheno_data)[1] <- "Sample_ID"
  colnames(pheno_data)[colnames(pheno_data) == "injection_order"] <- "Injection_order"
  
  mset <- construct_metabosets(exprs = exprs_data, pheno_data = pheno_data, feature_data = feature_data, group_col = "class", split_data = F)
  
  return(mset)
}
