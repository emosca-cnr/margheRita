#' Import a library from an msp file
#' @param msp_file msp file
#' @param name string for metabolite name among the spectraVariables()
#' @param rt string for precursor rt among the spectraVariables()
#' @param mz string for precursor mz among the spectraVariables()
#' @param CAS string for precursor CAS among the spectraVariables()
#' @param PubChemCID string for PubChemCID among the spectraVariables()
#' @param SMILES string for precursor SMILES among the spectraVariables()
#' @param accept_RI minimum relative intensity tof a peak to be kept in the library, default to 10
#' @importFrom Spectra Spectra spectraData peaksData
#' @importFrom MsBackendMsp MsBackendMsp
#' @export

library_from_msp <- function(msp_file=NULL, name="name", rt="rtime", mz="precursorMz", CAS=NULL, PubChemCID=NULL, SMILES=NULL, accept_RI=10){
	
	sps <- Spectra(msp_file,  source = MsBackendMsp(), backend = MsBackendMsp())
	
	ID_vec <- paste0("L", 1:length(sps))
	mR_lib <- list(
		lib_precursor=data.frame(ID=ID_vec, Name=spectraData(sps)[[name]], rt=spectraData(sps)[[rt]], mz=spectraData(sps)[[mz]], CAS=ifelse(is.null(CAS), NA, spectraData(sps)[[CAS]]), PubChemCID=ifelse(is.null(PubChemCID), NA, spectraData(sps)[[PubChemCID]]), SMILES=ifelse(is.null(SMILES), NA, spectraData(sps)[[SMILES]])),
		lib_peaks_data=data.frame(ID=ID_vec, Name=spectraData(sps)[[name]]),
		lib_peaks=setNames(lapply(peaksData(sps), function(pk) {
			colnames(pk) <- c("mz", "intensity")
			pk
		}), ID_vec),
		key_field="ID"
	)
	
	idx <- which(unlist(lapply(mR_lib$lib_peaks, function(x) nrow(x)>1)))
	mR_lib$lib_peaks <- mR_lib$lib_peaks[idx]
	
	mR_lib$lib_peaks <- lapply(mR_lib$lib_peaks, function(x){
		ans <- x
		ans[, 2] <- ans[, 2] / max(ans[, 2]) * 100
		ans <- ans[ans[, 2]>=accept_RI, , drop=F]
		return(ans)
	})
	
	mR_lib$lib_peaks_data <- mR_lib$lib_peaks_data[mR_lib$lib_peaks_data$ID %in% names(mR_lib$lib_peaks), ]
	
	return(mR_lib)
}