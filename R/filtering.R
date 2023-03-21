#' filtering
#' @param mRList mRList object
#' @param min_metab_in_sample min number of metabolites in a sample
#' @param min_sample_with_metab min number of samples in which a metabolite must appear
#' @param na_value value that indicate missing values
#' @param do_plot whether to plot or not
#' @param out_file output file name
#' @param lower_quality_mass_acc lower quality value
#' @param upper_quality_mass_acc upper quality value
#' @param color point color
#' @param ... further arguments to plot function
#' @return filtered mRList
#' @param seed set this seed with set.seed
#' @param a minimum increase factor
#' @param b maxmimum increase factor
#' @export
#' @importFrom stats density
#' @importFrom graphics polygon

#' @export
#' @return filtered mRList object

filtering <- function(mRList, seed=NULL, a=0.1, b=0.25, n=100, do_plot=TRUE, out_file="mz_quality.png", lower_quality_mass_acc=0.4, upper_quality_mass_acc=0.8, min_metab_in_sample=100, min_sample_with_metab=10, na_value="NA"){
 .filter_NA <- function(mRList, min_metab_in_sample=100, min_sample_with_metab=10, na_value="NA"){
    if(na_value != "NA"){
    cat("setting", na_value, "to NA\n")
    mRList$data[mRList$data == na_value] <- NA
  }

  idx_keep <- colSums(!is.na(mRList$data)) >= min_metab_in_sample
  cat("# Samples with >=", min_metab_in_sample, "metabolites", sum(idx_keep), "/", ncol(mRList$data), "\n")
  mRList$data <- mRList$data[, idx_keep]
  mRList$sample_ann <- mRList$sample_ann[idx_keep, ]

  idx_keep <- rowSums(!is.na(mRList$data)) >= min_sample_with_metab
  cat("# Features occurring in >= ", min_sample_with_metab, "samples", sum(idx_keep), "/", nrow(mRList$data), "\n")

  mRList$data <- mRList$data[idx_keep, ]
  mRList$metab_ann <- mRList$metab_ann[idx_keep, ]

  return(mRList)
}

.m_z_filtering <- function(mRList, do_plot=do_plot, out_file=out_file, lower_quality_mass_acc=lower_quality_mass_acc, upper_quality_mass_acc=upper_quality_mass_acc, color="black" , ...){

  mRList$metab_ann$quality <- quality <- mRList$metab_ann$mz %% 1

  all <- length(mRList$metab_ann$quality)
  idx_keep <- mRList$metab_ann$quality < lower_quality_mass_acc | mRList$metab_ann$quality > upper_quality_mass_acc

  cat("# Features with appropriate m/z values:", sum(idx_keep), "\n")
  cat("# Features without appropriate m/z values:", all-sum(idx_keep), "\n")

  mRList$data <- mRList$data[idx_keep, ]
  mRList$metab_ann <- mRList$metab_ann[idx_keep, ]

  mRList$QC <- mRList$QC[idx_keep, ]

  if(do_plot){

    grDevices::png(out_file, width = 200, height = 100, units="mm", res=300)
    par(mfrow=c(1, 2))

    plot(density(quality), xlab="m/z decimal part", main="Raw", ...)
    polygon(density(quality), col=color, border=color)


    plot(density(mRList$metab_ann$quality), xlab="m/z decimal part", main="Filtered", ...)
    polygon(density(mRList$metab_ann$quality), col=color, border=color)

    dev.off()

  }

  return(mRList)
}

.imputation <- function(mRList, seed=NULL, a=a, b=b, n=n) {
  
  
  #internal function
  imp_row <- function(row_min, n_na, a=a, b=b, n=100){
    
    ans <- seq(from = row_min * a, to = row_min * b, by = 1/n)
    ans <- sample(ans, n_na)
    return(ans)
    
  }
  
  if(!is.null(seed)){
    cat("Setting seed =", seed, ".\n")
    set.seed(seed)
  }
  
  na_val <- is.na(mRList$data)
  n <- min(n, max(rowSums(na_val))) #n is set possibly to the maximum number of na found in a row.
  
  #index of rows with NA
  idx_na_rows <- apply(na_val, 1, any)
  if(any(idx_na_rows)){
    
    #replacement
    mRList$data[idx_na_rows, ] <- apply(mRList$data[idx_na_rows, ], 1, function(x) replace(x = x, list = is.na(x), values = imp_row(min(x, na.rm = T), sum(is.na(x)), a=a, b=b, n=n)))
    
  }else{
    cat("Nothing to do.\n")
  }
  return(mRList)
}
 mRList <- filter_NA(mRList)
 mRList <- m_z_filtering(mRList)
 mRlist <- imputation(mRList)
 return(mRlist)
 }