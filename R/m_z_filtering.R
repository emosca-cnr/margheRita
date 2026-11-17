#' Filter features by the m/z decimal part
#' @param mRList mRList object
#' @param lower_quality_mass_acc lower quality value
#' @param upper_quality_mass_acc upper quality value
#' @return mRList object
#' @export

m_z_filtering <- function(mRList=NULL, lower_quality_mass_acc=0.4, upper_quality_mass_acc=0.8){
  
  mRList$metab_ann$quality <- quality <- mRList$metab_ann$mz %% 1
  
  all <- length(mRList$metab_ann$quality)
  idx_keep <- mRList$metab_ann$quality < lower_quality_mass_acc | mRList$metab_ann$quality > upper_quality_mass_acc
  
  cat("# Features with appropriate m/z values:", sum(idx_keep), "\n")
  cat("# Features without appropriate m/z values:", all-sum(idx_keep), "\n")
  
  mRList$data <- mRList$data[idx_keep, ]
  mRList$metab_ann <- mRList$metab_ann[idx_keep, ]
  
  mRList$QC <- mRList$QC[idx_keep, ]
  
  return(mRList)
}
