#' Transforms mRList into MSnSet object
#' @param mRList mRList object
#' @importFrom MSnbase MSnSet
#' @export
#' @return MSnSet object

as.MSnSet.mRList <- function(mRList){

  pheno_data <- rbind(mRList$sample_ann, mRList$QC_ann)
  pheno_data <- cbind(pheno_data, QC=pheno_data$class)

  msnset <- MSnbase::MSnSet(exprs = as.matrix(cbind(mRList$data, mRList$QC)), pData = pheno_data, fData = mRList$metab_ann)

  return(msnset)

}
