#' Convert sample intensity to relative intensity and filter
#'
#' @export
#' @param feature_spectra it is a list of sample ID which each contains list of its m/z and Intensity.
#' @param accept_RI numeric. the default value is 10. it is a maximum relative intensity that is kept in sample. since low intense peaks could be noise, it is filtering sample dataset by deleting the relative intensity lower then accept_RI.
#' @return list of spectra with RI

calc_RI <- function(feature_spectra=NULL, accept_RI = 10) {

  ans <- feature_spectra
  for (i in 1:length(ans)) {
    ans[[i]][, 2] <- ans[[i]][, 2] / max(ans[[i]][, 2]) * 100
    ans[[i]] <- ans[[i]][ans[[i]][, 2] >= accept_RI, , drop=FALSE]
  }
  
  return(ans)

}
