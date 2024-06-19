#' Filter m/z of metabolites
#' @description Function filters metabolites with 0.4>m/z<0.8 as it is not accurate.
#' @param mRList mRList object
#' @param lower_quality_mass_acc lower quality value
#' @param upper_quality_mass_acc upper quality value
#' @param color point color
#' @param dirout output directory
#' @param ... further arguments to function `plot()`
#' @export
#' @importFrom graphics hist
#' @return filteredd mRList object

m_z_filtering <- function(mRList=NULL, lower_quality_mass_acc=0.4, upper_quality_mass_acc=0.8, color="black" , dirout=NULL, ...){
  
  if (!is.null(dirout)) {
    dir.create(dirout, showWarnings = F)
  }else{
    dirout <- getwd()
  }
  
  
  mRList$metab_ann$quality <- quality <- mRList$metab_ann$mz %% 1
  
  all <- length(mRList$metab_ann$quality)
  idx_keep <- mRList$metab_ann$quality < lower_quality_mass_acc | mRList$metab_ann$quality > upper_quality_mass_acc
  
  cat("# Features with appropriate m/z values:", sum(idx_keep), "\n")
  cat("# Features without appropriate m/z values:", all-sum(idx_keep), "\n")
  
  mRList$data <- mRList$data[idx_keep, ]
  mRList$metab_ann <- mRList$metab_ann[idx_keep, ]
  
  mRList$QC <- mRList$QC[idx_keep, ]
  
  
  png(file.path(dirout, "hist.mz.png"), width = 200, height = 200, res=300, units="mm")
  par(mar=c(3, 3, .5, .5))
  par(mgp=c(1.7, .5, 0))
  #plot(density(quality), xlab="m/z decimal part", main="Raw", ...)
  #polygon(density(quality), col=color, border=color)
  hr <- hist(quality, breaks=20, plot=F, ...)
  hist(quality, breaks=20, xlab="m/z decimal part", main="", freq=F, ylim=c(0, max(hr$density)+1), col="dodgerblue", ...)
  dev.off()
   
  png(file.path(dirout, "hist.mz.filtered.png"), width = 200, height = 200, res=300, units="mm")
  par(mar=c(3, 3, .5, .5))
  par(mgp=c(1.7, .5, 0))
  #plot(density(mRList$metab_ann$quality), xlab="m/z decimal part", main="Filtered", ...)
  #polygon(density(mRList$metab_ann$quality), col=color, border=color)
  hr <- hist(mRList$metab_ann$quality, breaks=20, plot=F)
  hist(mRList$metab_ann$quality, breaks=20, xlab="m/z decimal part", main="", freq=F, ylim=c(0, max(hr$density)+1), col="dodgerblue", ...)
  dev.off()
  
  return(mRList)
}
