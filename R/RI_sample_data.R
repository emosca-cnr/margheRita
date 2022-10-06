#' Convert sample intensity to relative intensity and filter
#'
#' @export
#' @param feature_spectra it is a list of sample ID which each contains list of its m/z and Intensity.
#' @param accept_RI numeric parameter. the default value is 10. it is a maximum relative intensity that is kept in sample. since low intense peaks could be noise, it is filtering sample dataset by deleting the relative intensity lower then accept_RI.
#'
RI_sample_data = function(feature_spectra, accept_RI = 10) {

  #calculating relative intensity
  RI_sample1 = lapply(1 :length(feature_spectra), function(x) feature_spectra[[x]][,2] / max(feature_spectra[[x]][,2]) * 100 )

  #establishing new data set for sample spectra
  names(RI_sample1) = names(feature_spectra)
  RI_sample  = feature_spectra

  for (n in 1:length(RI_sample)) {
    RI_sample[[n]][,2] = RI_sample1[[n]]
  }

  # selected Relative intensity > 10 with correlated mass in sample
  for (m in 1: length(RI_sample)) {
    RI_sample[[m]] = RI_sample[[m]][RI_sample[[m]][,2]>= accept_RI , , drop=FALSE]
  }

  return(RI_sample)

}
