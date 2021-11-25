#' establishing new sample data by calculating relative intensity
#'
#'
#' @export
#'

RI_sample_data = function(feature_data_spectra ) {

# selected Relative intensity > 15 with correlated mass in sample
RI_sample1 = lapply(1 : length(feature_data_spectra), function(x) feature_data_spectra[[x]]$i / max(feature_data_spectra[[x]]$i) * 100 )
names(RI_sample1) = names(feature_data_spectra)
RI_sample  = feature_data_spectra
for (n in 1:length(RI_sample)) { RI_sample[[n]][,2] = RI_sample1[[n]]}
for (m in 1: length(RI_sample)) {RI_sample[[m]] = RI_sample[[m]][RI_sample[[m]][,2]>=15, ]}

return(RI_sample)

}
