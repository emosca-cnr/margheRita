#' Convert spectra from character vector to a list of data frames
#' @param spectra a character vector where each element is a spectra, encoded as mz1:I1 mz2:I2, for example: "45.9732:9894 46.97655:21358..."
#' @export
#' @return a list of spectra, where each spectra is a data.frame with two columns (mz and intensity)

get_spectra_list_from_vector <- function(spectra=NULL){

	ref_lib_spectra <- sapply(spectra, function(x) strsplit(x, " "))
	ref_lib_spectra <- lapply(ref_lib_spectra, function(x) lapply(strsplit(x, ":"), as.numeric))
	ref_lib_spectra <- lapply(ref_lib_spectra, function(x) do.call(rbind, x))
	names(ref_lib_spectra) <- names(spectra)

	return(ref_lib_spectra)


}
