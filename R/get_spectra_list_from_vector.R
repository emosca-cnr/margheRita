#'get_spectra_list_from_vector
#' @export

get_spectra_list_from_vector <- function(spectra){

	ref_lib_spectra <- sapply(spectra, function(x) strsplit(x, " "))
	ref_lib_spectra <- lapply(ref_lib_spectra, function(x) lapply(strsplit(x, ":"), as.numeric))
	ref_lib_spectra <- lapply(ref_lib_spectra, function(x) do.call(rbind, x))
	names(ref_lib_spectra) <- names(spectra)

	return(ref_lib_spectra)


}
