#' m_z_filtering
#' @param mRList
#' @export

m_z_filtering <- function(mRList, do_plot=TRUE, lower_quality_mass_acc=0.4, upper_quality_mass_acc=0.8, color="black" ,...){

  mRList$metab_ann$quality <- mRList$metab_ann$mz %% 1
  all <- length(mRList$metab_ann$quality)
  idx_keep <- mRList$metab_ann$quality < lower_quality_mass_acc | mRList$metab_ann$quality > upper_quality_mass_acc

  cat("Metabolites with appropriate m/z values\n")
  cat(as.data.frame(table(idx_keep))[1, 2])
  cat("\n")
  cat("Metabolites witout appropriate m/z values\n")
  cat(all-(as.data.frame(table(idx_keep)))[1, 2])

  mRList$data <- mRList$data[idx_keep, ]
  mRList$metab_ann <- mRList$metab_ann[idx_keep, ]

  mRList$QC <- mRList$QC[idx_keep, ]

   if(do_plot){
    plot(density(mRList$metab_ann$quality), xlab="m/z Quality", main="Distrubution of m/z Quality", ...)+
    polygon(density(mRList$metab_ann$quality), col=color, border=color)
  }
  return(mRList)
}
