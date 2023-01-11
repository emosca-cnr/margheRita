#' Library selection
#' Library of margheRita contains three list:
#' - precursors library: it is a list of metabolites which each contains information of retention time and mz with their specific IDs, CAS numbers and names.
#' - lib_peaks_data:it is a list of metabolites which each contains information of collision energy with their specific IDs, CAS numbers and names.
#' - lib_peaks: it is a list of metabolites with same IDs as lib_peaks_data which each contains a list of a peaks with mz and relative intensity.
#'
#' @param column column type, only if source=="margheRita". Possible values are: HILIC, LipC8, pZIC, RPLong, RPShort. Based on the type of the column, the retention time of same metabolite could be different.
#' @param mode mode could be set in positive or negative state. positive mode select positive collision energy and mz in positive mode.
#' @param accept_RI numeric parameter. the default value is 10. it is a maximum relative intensity that we keep in library. since low intense peaks could be noise, it is filtering the library by deleting the relative intensity lower then accept_RI.
#'
#' @export
#' @importFrom utils data
select_library <- function(source=c("margheRita", "MS-DIAL"), column=NULL, mode=c("POS", "NEG"), accept_RI=10){
  
  mRlib_peaks_list <- mRlib_peaks_df <- mRlib_precursors <- NULL #to please the check
  
  source <- match.arg(source, c("margheRita", "MS-DIAL"))
  column <- match.arg(column, c("HILIC", "LipC8", "pZIC", "RPLong", "RPShort"))
  mode <- match.arg(mode, c("POS", "NEG"))
  
  cat("Library source:", source, "\n")
  if(source == "margheRita"){
    cat("Column:", column, "\n")
  }
  cat("mode: ", mode, "\n")
  cat("minimum RI: ", accept_RI, "\n")
  
  if(source=="margheRita"){
    
    data("mRlib_peaks_list", envir=environment())
    data("mRlib_peaks_df", envir=environment())
    data("mRlib_precursors", envir=environment())
    
    lib_precursor <- unique(mRlib_precursors[mRlib_precursors$COL==column, c("ID", "CAS", "Name", "rt", mode, "PubChemCID")])
    rownames(lib_precursor) <- lib_precursor$ID
    colnames(lib_precursor)[5] <- "mz"
    
    lib_peaks_data <- unique(mRlib_peaks_df[, c("ID", "CAS", "Collision_energy", "Name")])
    lib_peaks_data$Collision_energy[lib_peaks_data$Collision_energy < 0] <- "NEG"
    lib_peaks_data$Collision_energy[lib_peaks_data$Collision_energy != "NEG"] <- "POS"
    lib_peaks_data <- lib_peaks_data[lib_peaks_data$Collision_energy == mode, ]
    
    lib_peaks <- mRlib_peaks_list
    
    key_field <- "CAS"
  }
  
  if(source=="MS-DIAL"){
    
    input_data_file <- system.file("extdata", paste0("MSDIAL_precursors_", mode, ".rds"), package = "margheRita")
    lib_precursor <- readRDS(input_data_file)
    rownames(lib_precursor) <- lib_precursor$ID
    
    input_data_file <- system.file("extdata", paste0("MSDIAL_peaks_df_", mode, ".rds"), package = "margheRita")
    lib_peaks_data <- readRDS(input_data_file)
    
    input_data_file <- system.file("extdata", paste0("MSDIAL_peaks_list_", mode, ".rds"), package = "margheRita")
    lib_peaks <- readRDS(input_data_file)
    
    key_field <- "ID"
  }
 
  #filter peaks by RI
  lib_peaks <- lapply(lib_peaks, function(x) x$peaks[x$peaks[, 3] >= accept_RI, c(1, 3), drop=FALSE])
 
  #match the elements of the objects related to 
  lib_peaks <- lib_peaks[match(lib_peaks_data$ID, names(lib_peaks))]
  
  cat("# metabolites: ", nrow(lib_precursor), "\n")

  return(list(lib_precursor=lib_precursor, lib_peaks=lib_peaks, lib_peaks_data=lib_peaks_data, key_field=key_field))
  
}
