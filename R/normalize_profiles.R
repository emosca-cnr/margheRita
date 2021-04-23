#' Normalize profiles
#' @param m_list$df features-by-samples matrix of metabolite levels;
#' @param method "log": log2(x+1); "reference": divide each sample by its reference; "pqn": probailistic quotient normalization;
#' @param reference reference profile;
#' @export

normalize_profiles <- function(m_list, method=c("log", "reference", "pqn")){

  if(method == "log"){
    m_list$df <- log2(m_list$df+1)
  }

  if(method == "reference"){
    m_list$df <- t(t(m_list$df) / m_list$metab_ann$reference)
  }

  if(method == "pqn"){
    m_list$df <- apply(m_list$df, 2, function(x) pqn(x, xref = m_list$metab_ann$reference))
  }

  return(m_list)

}
