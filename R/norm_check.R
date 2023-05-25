#' Test the normality of each metabolite by Shapiro-Wilk's test
#'
#' Function calculates p-values for each metabolites by Shapiro-Wilk's test to test normal distribution. If p-value>0.05 metabolite is normally distributed, p-value<0.05 metabolite is not normally distributed
#' @param mRList mRList object
#' @return mRList object with mRist$Shapiro element that reports p-value calculated using Shapiro-Wilk's test
#' @export
#' @importFrom stats shapiro.test



norm_check <- function(mRList=NULL) {

  Shapirotest <- apply(mRList$data, 1, function(x) shapiro.test(as.numeric(x))$"p.value"[1])
  Shapirotest<-as.data.frame(Shapirotest)
  rownames(Shapirotest) <- rownames(mRList$data)
  mRList$shapirotest <- Shapirotest

  return(mRList)
}

