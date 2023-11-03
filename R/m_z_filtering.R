#' Filter m/z of metabolites
#' @description Function filters metabolites with 0.4>m/z<0.8 as it is not accurate.
#' @param mRList mRList object
#' @param do_plot whether to plot or not
#' @param out_file output file name
#' @param lower_quality_mass_acc lower quality value
#' @param upper_quality_mass_acc upper quality value
#' @param color point color
#' @param ... further arguments to plot function
#' @return filtered mRList
#' @export
#' @importFrom stats density
#' @importFrom graphics polygon
#' @importFrom grDevices png
#' @return mRList object filtered by m/z not accurate
#' @examples
#' ##library(dataset.margheRita)
#' ##dataset(norm_pos)
#'mRList<-m_z_filtering(mRList,do_plot=TRUE, out_file="mz_quality.png", lower_quality_mass_acc=0.4, upper_quality_mass_acc=0.8, color="black")
#'

m_z_filtering <- function(mRList, do_plot=TRUE, lower_quality_mass_acc=0.4, upper_quality_mass_acc=0.8, color="black" , ...){
  
  mRList$metab_ann$quality <- quality <- mRList$metab_ann$mz %% 1
  
  all <- length(mRList$metab_ann$quality)
  idx_keep <- mRList$metab_ann$quality < lower_quality_mass_acc | mRList$metab_ann$quality > upper_quality_mass_acc
  
  cat("# Features with appropriate m/z values:", sum(idx_keep), "\n")
  cat("# Features without appropriate m/z values:", all-sum(idx_keep), "\n")
  
  mRList$data <- mRList$data[idx_keep, ]
  mRList$metab_ann <- mRList$metab_ann[idx_keep, ]
  
  mRList$QC <- mRList$QC[idx_keep, ]
  
  
  plot(density(quality), xlab="m/z decimal part", main="Raw", ...)
  polygon(density(quality), col=color, border=color)
  
  plot(density(mRList$metab_ann$quality), xlab="m/z decimal part", main="Filtered", ...)
  polygon(density(mRList$metab_ann$quality), col=color, border=color)
  
  
  
  return(mRList)
}
