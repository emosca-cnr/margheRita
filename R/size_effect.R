#' size effect for pqn
#' @param x metabolomics profile
#' @param xref reference profile
#' @param robust TRUE/FALSE
#' @value list of:
#'    - se = size effects;
#'    - med = average size effect.
#' @author Ettore Mosca

size_effect <- function(x, xref, robust=T){

  #full size effects
  se <- x / xref

  #median ration between elements of x and xref
  if(robust){
  	se_med <- median(se[which(se>0)], na.rm = T)
  }else{
  	se_med <- mean(se[which(se>0)], na.rm = T)
  }

  return(list(se=se, med=se_med))

}
