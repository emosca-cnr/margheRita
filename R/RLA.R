#' Relative Log Abudance
#' @param x metabolites-by-samples matrix of metabolite levels
#' @export
#' @author Ettore Mosca (CNR-ITB)
#' @importFrom stats median

RLA <- function(x, robust=TRUE){

	ans <- log(x+1)
	if(robust){
		ans <- t(apply(ans, 1, function(y) y - median(y)))
	}else{
		ans <- t(apply(ans, 1, function(y) y - mean(y)))
	}

	return(ans)

}