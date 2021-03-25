#' size effect for pqn
#' @param x metabolomics profile
#' @param xref reference profile
#' @param robust TRUE/FALSE
#' @importFrom stats median
#' @export

#' @return list of:
#'    - se = size effects;
#'    - avg = average size effect.
#' @author Ettore Mosca (CNR-ITB)

size_effect <- function(x, xref, robust=T){

  #full size effects
  se <- x / xref

  #median ration between elements of x and xref
  if(robust){
  	se_avg <- median(se[which(se>0)], na.rm = T)
  }else{
  	se_avg <- mean(se[which(se>0)], na.rm = T)
  }

  return(list(se=se, avg=se_avg))

}
