#' m_z_filtering
#' @param mRList mRList
#' @export

m_z_filtering <- function(mRList, do_plot=TRUE, out_file="mz_quality.png", lower_quality_mass_acc=0.4, upper_quality_mass_acc=0.8, color="black" ,...){
  
  mRList$metab_ann$quality <- mRList$metab_ann$mz %% 1
  all <- length(mRList$metab_ann$quality)
  idx_keep <- mRList$metab_ann$quality < lower_quality_mass_acc | mRList$metab_ann$quality > upper_quality_mass_acc
  
  cat("# Metabolites with appropriate m/z values:", sum(idx_keep), "\n")
  cat("# Metabolites witout appropriate m/z values:", all-sum(idx_keep), "\n")
  
  mRList$data <- mRList$data[idx_keep, ]
  mRList$metab_ann <- mRList$metab_ann[idx_keep, ]
  
  mRList$QC <- mRList$QC[idx_keep, ]
  
  if(do_plot){
    
    grDevices::png(out_file, width = 180, height = 180, units="mm", res=300)
    
    plot(density(mRList$metab_ann$quality), xlab="m/z decimal part", main="Distrubution of m/z Quality", ...)
    polygon(density(mRList$metab_ann$quality), col=color, border=color)
    
    dev.off()
    
  }
  
  return(mRList)
}
