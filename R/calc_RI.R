#' Convert sample intensity to relative intensity and filter
#'
#' @param feature_spectra list of spectra; each spectrum is a matrix with two columns in the mandatory order: m/z, intensity.
#' @param accept_RI numeric, default to 10. The minimum value of relative intensity that will be kept; the peaks with R less than accept_RI are removed.
#' @return list of spectra where intensity is replaced by RI

calc_RI <- function(feature_spectra=NULL, accept_RI = 10) {

  ans <- feature_spectra
  for (i in 1:length(ans)) {
    ans[[i]][, 2] <- ans[[i]][, 2] / max(ans[[i]][, 2]) * 100
    ans[[i]] <- ans[[i]][ans[[i]][, 2] >= accept_RI, , drop=FALSE]
  }
  
  return(ans)

}
