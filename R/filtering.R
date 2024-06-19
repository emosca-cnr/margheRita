#' Filtering of features
#' @description Filtering of features by missing values, m/z values and imputation of missing values
#' @param mRList mRList object
#' @param min_metab_in_sample min number of metabolites in a sample
#' @param min_sample_with_metab min number of samples in which a metabolite must appear
#' @param na_value value that indicate missing values
#' @param lower_quality_mass_acc lower quality value
#' @param upper_quality_mass_acc upper quality value
#' @param seed set this seed with set.seed
#' @param a minimum increase factor
#' @param b maxmimum increase factor
#' @param dirout output directory
#' @export
#' @return filtered mRList object

filtering <- function(mRList=NULL, seed=NULL, a=0.1, b=0.25, lower_quality_mass_acc=0.4, upper_quality_mass_acc=0.8, min_metab_in_sample=100, min_sample_with_metab=3, na_value="NA", dirout=NULL){
  
  mRList <- filter_NA(mRList = mRList, min_metab_in_sample = min_metab_in_sample, min_sample_with_metab = min_sample_with_metab, na_value = na_value)
  
  mRList <- m_z_filtering(mRList = mRList, lower_quality_mass_acc = lower_quality_mass_acc, upper_quality_mass_acc = upper_quality_mass_acc, dirout=dirout)
  
  mRList <- imputation(mRList = mRList, seed = seed, a = a, b=b)
  
  return(mRList)
}