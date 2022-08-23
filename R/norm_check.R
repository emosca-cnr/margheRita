#' Test the normality of each metabolite by Shapiro-Wilk's test
#' @description Function calculates p-values for each metabolites by Shapiro-Wilk's test to test normal distribution.
#' @description If p-value>0.05 metabolite is normally distributed, p-value<0.05 metabolite is not normally distributed
#' @param mRList mRList object
#' @return mRist$Shapiro that reports p-value calculated using Shapiro-Wilk's test
#' @export



norm_check <- function(mRList) {

  Shapirotest <- apply(mRList$data,1,function(x) shapiro.test(as.numeric(x)))
  mRList$shapirotest<-Shapirotest

  return(mRList)
}

