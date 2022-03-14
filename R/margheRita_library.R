#' Library of margheRita
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
#'@return
#'@export
#'
#'@examples
#'
margheRita_library <- function(column=c("HILIC", "LipC8", "pZIC", "RPLong", "RPShort"), mode=c("POS", "NEG"), accept_RI = 10){
  
  ## library
  data("mRlib_peaks_df", envir=environment())
  data("mRlib_precursors", envir=environment())
  
  
  lib_precursors <- mRlib_precursors[mRlib_precursors$COL==column, c( "ID", "CAS", "Name", "rt", mode)]
  rownames(lib_precursors) <- lib_precursors$ID
  colnames(lib_precursors)[5] <- "mz"
  
  lib_peaks_data <- unique(mRlib_peaks_df[ , c("ID", "CAS", "Collision_energy", "Name", "Peaks")])
  lib_peaks <- sapply(lib_peaks_data$Peaks, function(x) strsplit(x, " "))
  lib_peaks <- lapply(lib_peaks, function(x) lapply(strsplit(x, ":"), as.numeric))
  lib_peaks <- lapply(lib_peaks, function(x) do.call(rbind, x))
  names(lib_peaks) <- lib_peaks_data$ID
  
  for(e in 1:length(lib_peaks)) {
    lib_peaks[[e]] <- lib_peaks[[e]][lib_peaks[[e]][,2] > accept_RI, , drop=FALSE]
  }
  
  lib_peaks_data$Peaks <- NULL
  
  lib_peaks_data$Collision_energy[lib_peaks_data$Collision_energy < 0] <- "NEG"
  lib_peaks_data$Collision_energy[lib_peaks_data$Collision_energy != "NEG"] <- "POS"
  lib_peaks_data <- lib_peaks_data[lib_peaks_data$Collision_energy == mode, ]
  
  
  return(list(lib_precursors = lib_precursors, lib_peaks=lib_peaks, lib_peaks_data=lib_peaks_data))
  
}