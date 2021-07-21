#' Shapiro-Wilk's test is used for testing the normal distribution of metabolites lo faccio sul subsetting Samples

#' @param m_list margheRita m_list
#' @export



norm_check <- function(m_list) {

  Shapirotest <- apply(m_list$data,1,function(x) shapiro.test(as.numeric(x)))
  m_list$shapirotest<-Shapirotest

  return(m_list)
}

