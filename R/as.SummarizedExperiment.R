#' from mRList to SummarizedExperiment 
#' @param mRList mRList object
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @return object of class SummarizedExperiment

as.SummarizedExperiment <- function(mRList=NULL){
  
  ### 1.1 #### create a Summarized experiment
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
  
  se <-  SummarizedExperiment(
    assays = exprs_data,
    rowData = feature_data,
    colData = pheno_data
  )
  
  
  return(se)
}
