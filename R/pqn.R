#' PQN normalization
#' @param x metabolomics profile
#' @param xref reference profile
#' @return list of:
#'    - y = normalized profile;
#'    - se = size effects;
#'    - se_med = average size effect.
#' @export
#' @author Ettore Mosca


pqn <- function(x, xref){

	#size effetcs:
	#se$se: full size effect matrix
	#se$med: #medians of size effects (calculated for positive elements)
	se <- size_effect(x, xref)

  #pqn
  out <- x / se$avg

  return(list(y=out, se=se$se, se_avg=se$avg))

}
