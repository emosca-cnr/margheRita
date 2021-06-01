#' Normalize profiles
#' @param m_list$data features-by-samples matrix of metabolite levels;
#' @param method "log": log2(x+1); "reference": divide each sample by its reference; "pqn": probailistic quotient normalization;
#' @param reference reference profile;
#' @export

normalize_profiles <- function(m_list, method=c("log", "reference", "pqn")){

  if(method == "log"){
    m_list$data <- log2(m_list$data+1)
  }

  if(method == "reference"){
    m_list$data <- t(t(m_list$data) / m_list$metab_ann$reference)
  }

  if(method == "pqn"){
    m_list$data <- apply(m_list$data, 2, function(x) pqn(x, xref = m_list$metab_ann$reference)$y)
  }

  return(m_list)

}
