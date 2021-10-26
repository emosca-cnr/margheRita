#' Shapiro-Wilk's test is used for testing the normal distribution of metabolites
#' Calculate p-value, if p-value>0.05 metabolite is normally distributed, p-value<0.05 metabolite is not normally distributed
#' @param m_list
#' @return  m_list Shapiro with p-value calculated usign Shapiro-Wilk's test
#' @export



norm_check <- function(m_list) {

  Shapirotest <- apply(m_list$data,1,function(x) shapiro.test(as.numeric(x)))
  m_list$shapirotest<-Shapirotest

  return(m_list)
}

