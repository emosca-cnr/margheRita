#' Transforms mRList into PomaSummarizedExperiment
#' @param mRList mRList object
#' @export
#' @return a list that cotains the arguments for notame::construct_metabosets function

as.PomaSummarizedExperiment <- function(mRList=NULL){
  
  pheno_data <- mRList$sample_ann
  exprs_data <- mRList$data
  
  if(length(mRList$QC_ann) > 0){
    pheno_data <- rbind(mRList$sample_ann, mRList$QC_ann)
    pheno_data <- cbind(pheno_data, QC=pheno_data$class)
    exprs_data <- cbind(exprs_data, mRList$QC)
  }
  colnames(pheno_data)[1] <- "Sample_ID"
  colnames(pheno_data)[colnames(pheno_data) == "injection_order"] <- "Injection_order"
  
  pse <- POMA::PomaSummarizedExperiment(target = pheno_data, features = exprs_data)
  
  return(pse)
}
