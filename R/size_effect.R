#' Size effect for pqn
#' @param x metabolomics profile
#' @param xref reference profile
#' @param robust TRUE/FALSE
#' @importFrom stats median
#' @return list of:
#'    - se = size effects;
#'    - avg = average size effect.
#' @export

size_effect <- function(x=NULL, xref=NULL, robust=T){

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
