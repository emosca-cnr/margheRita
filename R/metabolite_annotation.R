#' Metabolite Annotation
#' Run the metabolite identification pipeline
#' @param mRList mRList object
#' @param library_list library data obtained from margherita_library() function
#' @param rt_err threshold over the retention time error
#' @param accept_flag A number with default value of 5. PPM errors < accept_flag will be tagged as "super", while those > accept_flag and < suffer_flag will be tagged as "acceptable"
#' @param suffer_flag A number with default value of 10. PM errors above this value and < unaccept_flag will be tagged as "suffer"
#' @param unaccept_flag A number with default value of 15. The maximum PPM error must be less than this value. and those above this number will be eliminated.
#' @param min_RI numeric parameter. the default value is 10. it is a maximum relative intensity that is kept in sample. since low intense 
#' peaks could be noise, it is filtering sample dataset by deleting the relative intensity lower then accept_RI.
#' @param ppm_err A number with default value of 10. The maximum PPM error must be less than this value. and those above this number will be eliminated.
#' @param mode mode could be set in positive or negative state. positive mode select positive collision energy and mz in positive mode.
#' @param RI_err maximum absolute RI difference between MS/MS peaks of sample and library
#' @param RI_err_type type of RI error calculation.

#' @return data.frame with matches between features and library metabolites

#' @export
#' @importFrom stats setNames

metabolite_annotation <- function(mRList = NULL, library_list = NULL, rt_err=1, unaccept_flag=20, accept_flag=5, suffer_flag=10, min_RI = 10, ppm_err = 20, mode=NULL, RI_err=20, RI_err_type="rel"){
  
  #extract feature information from metabolite metadata
  out_levels <- unique(mRList_raw$metab_ann[, c("Feature_ID", "rt", "mz", "MS_MS_spectrum")])
  
 
  ### check retention time and mass
  cat("Checking RT... ")
  RT <-  check_RT(feature_data = out_levels, reference = library_list$lib_precursor, rt_err_thr=rt_err)
  cat("done\n")
  
  cat("Checking precursor mz... ")
  mass <- check_mass(feature_data = out_levels, reference = library_list$lib_precursor, unaccept_flag=unaccept_flag, accept_flag=accept_flag, suffer_flag=suffer_flag)
  cat("done\n")
  
  cat("Integrating RT and mz results... ")
  RT_mass <- check_RT_mass(RT = RT, mass = mass)
  cat("done\n")
  
  cat("Assembling output for precursors... ")
  ### features with appropriate mass
  mass_ok <- stack(lapply(mass, function(x) x$Feature_ID))
  mass_ok <- merge(library_list$lib_precursor, mass_ok, by.x="ID", by.y="ind", all.y=T)
  colnames(mass_ok) <- replace(colnames(mass_ok), list = colnames(mass_ok)=="values", "Feauture_ID")
  
  out_levels <- merge(out_levels, data.frame(ID=rep(names(mass), unlist(lapply(mass, nrow))), do.call(rbind, lapply(mass, function(x) x[, colnames(x) != "mz"])), stringsAsFactors = F), by="Feature_ID", all=T)
  out_levels <- merge(library_list$lib_precursor, out_levels, by="ID", all = T, suffixes = c("_lib", ""))
  
  ### features with appropriate mass and RT
  RT_mass_ok <- stack(lapply(RT_mass, function(x) x$Feature_ID))
  RT_mass_ok <- merge(library_list$lib_precursor, RT_mass_ok, by.x="ID", by.y="ind", all.y=T)

  out_levels <- merge(out_levels, data.frame(ID=rep(names(RT_mass), unlist(lapply(RT_mass, nrow))), do.call(rbind, lapply(RT_mass, function(x) x[, c("Feature_ID", "rt", "RT_err", "RT_flag")])), stringsAsFactors = F), by=c("ID", "Feature_ID", "rt"), all=T)
  cat("done\n")
  
  
  ###MS/MS matching - at least mass
  cat("MS/MS analysis... ")
  feature_spectra_list <- get_spectra_list_from_vector(spectra = setNames(out_levels$MS_MS_spectrum[!is.na(out_levels$MS_MS_spectrum)], out_levels$Feature_ID[!is.na(out_levels$MS_MS_spectrum)]))
  RI_sample <- RI_sample_data(feature_spectra=feature_spectra_list, accept_RI = min_RI)
  
  intense_peak <- peak_matching(RT_mass = mass, RI_lib = library_list$lib_peaks, reference = library_list$lib_precursor, lib_peaks_data=library_list$lib_peaks_data, mode = mode, RI_sample = RI_sample, ppm_err = ppm_err, intensity = RI_err, RI_diff_type = RI_err_type)
  cat("done\n")
  
  
  cat("Assembling final output...  ")
  
  intense_peak_unq <- unique(intense_peak$matched_peaks[, c("ID", "Feature_ID")])
  intense_peak_unq <- merge(library_list$lib_precursor, intense_peak_unq, by="ID", all.y=T)
  

  out_levels <- merge(out_levels, intense_peak$matched_peaks[, c("ID", "Feature_ID", "ID_peaks", "peaks_found_ppm_RI", "precursor_in_MSMS")], by=c("ID", "Feature_ID"), all=T)
  
  out_levels$Level <- ""
  out_levels$Level_note <- ""
  
  idx <- which(out_levels$mass_flag)
  out_levels$Level[idx] <- 3
  out_levels$Level_note[idx] <- "mz"
  
  idx <- which(out_levels$mass_flag & out_levels$RT_flag)
  out_levels$Level[idx] <- 3
  out_levels$Level_note[idx] <- "mz, rt"
  
  idx <- which(!is.na(out_levels$peaks_found_ppm_RI))
  out_levels$Level[idx] <- 2
  
  idx <- which(!is.na(out_levels$peaks_found_ppm_RI) & out_levels$mass_flag & out_levels$RT_flag)
  out_levels$Level[idx] <- 1
  
  
  ### summary
  out_levels <- out_levels[out_levels$Level != "", ]
  out_levels <- out_levels[, c("Feature_ID", "rt", "mz", "MS_MS_spectrum", "RT_err", "RT_flag", "ppm_error", "mass_flag", "mass_status", "ID_peaks", "peaks_found_ppm_RI", "precursor_in_MSMS", "ID", "Name", "rt_lib", "mz_lib", "Level", "Level_note", "CAS", "PubChemCID")]
  
  cat("done\n")
  
  return(out_levels)
  
  
 }
