#' establishing new sample data by calculating relative intensity
#'
#'
#' @export
#'

RI_sample_data = function(feature_spectra, acceptable_RI = acceptable_RI) {

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
  RI_sample[[m]] = RI_sample[[m]][RI_sample[[m]][,2]>= acceptable_RI , , drop=FALSE]
  }

return(RI_sample)

}
