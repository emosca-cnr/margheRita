#' Library selection
#' @export
#'
#' Library of margheRita contains three list:
#' - precursors library: it is a list of metabolites which each contains information of retention time and mz with their specific IDs, CAS numbers and names.
#' - lib_peaks_data:it is a list of metabolites which each contains information of collision energy with their specific IDs, CAS numbers and names.
#' - lib_peaks: it is a list of metabolites with same IDs as lib_peaks_data which each contains a list of a peaks with mz and relative intensity.
#'
#'
#'@param column it is specify the type of the column which could be; HILIC, LipC8, pZIC, RPLong, RPShort. based on the type of the column, the retention time of same metabolite could be different.
#'@param mode mode could be set in positive or negative state. positive mode select positive collision energy and mz in positive mode.
#'@param accept_RI numeric parameter. the default value is 10. it is a maximum relative intensity that we keep in library. since low intense peaks could be noise, it is filtering the library by deleting the relative intensity lower then accept_RI.
#'
select_library <- function(column=c("HILIC", "LipC8", "pZIC", "RPLong", "RPShort"), mode=c("POS", "NEG"), accept_RI=10){

  data("mRlib_peaks_list", envir=environment())
  data("mRlib_peaks_df", envir=environment())
  data("mRlib_precursors", envir=environment())

  lib_precursor <- unique(mRlib_precursors[mRlib_precursors$COL==column, c("ID", "CAS", "Name", "rt", mode)])
  rownames(lib_precursor) <- lib_precursor$ID
  colnames(lib_precursor)[5] <- "mz"

  #we need to use ID, because we have duplicated names with different peaks
  lib_peaks <- lapply(mRlib_peaks_list, function(x) x$peaks[x$peaks[, 3] > accept_RI, c(1, 3), drop=FALSE])

  lib_peaks_data <- unique(mRlib_peaks_df[, c("ID", "CAS", "Collision_energy", "Name")])
  lib_peaks_data$Collision_energy[lib_peaks_data$Collision_energy < 0] <- "NEG"
  lib_peaks_data$Collision_energy[lib_peaks_data$Collision_energy != "NEG"] <- "POS"
  lib_peaks_data <- lib_peaks_data[lib_peaks_data$Collision_energy == mode, ]

  lib_peaks <- lib_peaks[match(lib_peaks_data$ID, names(lib_peaks))]

  return(list(lib_precursor=lib_precursor, lib_peaks=lib_peaks, lib_peaks_data=lib_peaks_data))

}
