#' as.MSnSet.mRList
#'
#'
#' @importFrom MSnbase MSnSet
#' @export

as.MSnSet.mRList <- function(mRlist){

  pheno_data <- rbind(mRList$sample_ann, mRList$QC_ann)
  pheno_data <- cbind(pheno_data, QC=pheno_data$class)

  msnset <- MSnbase::MSnSet(exprs = as.matrix(cbind(mRlist$data, mRlist$QC)), pData = pheno_data, fData = mRList$metab_ann)

  return(msnset)

}
