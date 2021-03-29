#' Normalize profiles
#' @param X features-by-samples matrix of metabolite levels;
#' @param method "log": log-normalization; "reference": divide each sample by its reference; "pqn": probailistic quotient normalization;
#' @param reference reference profile;
#' @export

normalize_profiles <- function(X, method=c("log", "reference", "pqn"), reference=NULL){

  if(method == "log"){
    ans <- log2(X)
  }

  if(method == "reference"){
    ans <- t(t(X) / reference)
  }

  if(method == "pqn"){
    ans <- apply(X, 2, function(x) pqn(x, xref = reference))
  }

  return(ans)

}
