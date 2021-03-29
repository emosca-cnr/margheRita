#' Normalize profiles
#' @param X features-by-samples matrix of metabolite levels
#' @param method log log-normalization; reference: divide each sample by its reference; pqn probailistic quotient normalization #' @param pqn_reference reference profile for PQN
#' @export

normalize_profiles <- function(X, method=c("log", "reference", "pqn"), pqn_reference=NULL){


  if(method == "log"){
    ans <- log2(X)
  }

  if(method == "reference"){
    ans <- t(t(X) / totproteins$TOT_PROTEINS)
  }

  if(method == "pqn"){
    ans <- apply(X, 2, function(x) pqn(x, xref = pqn_reference))
  }

  return(ans)

}
