#' Metabolite Annotation
#' Run the metabolite identification pipeline
#' @param mRList mRList object
#' @param library_list library data obtained from margherita_library() function
#' @param feature_spectra MS_MS spectra (optional). By default spectra are taken from mRList object
#' @param rt_err_thr threshold over the retention time error
#' @param accept_flag A number with default value of 5. PPM errors < accept_flag will be tagged as "super", while those > accept_flag and < suffer_flag will be tagged as "acceptable"
#' @param suffer_flag A number with default value of 10. PM errors above this value and < unaccept_flag will be tagged as "suffer"
#' @param unaccept_flag A number with default value of 15. The maximum PPM error must be less than this value. and those above this number will be eliminated.
#' @param acceptable_RI numeric parameter. the default value is 10. it is a maximum relative intensity that is kept in sample. since low intense 
#' peaks could be noise, it is filtering sample dataset by deleting the relative intensity lower then accept_RI.
#' @param acceptable_PPM_err A number with default value of 10. The maximum PPM error must be less than this value. and those above this number will be eliminated.
#' @param mode mode could be set in positive or negative state. positive mode select positive collision energy and mz in positive mode.
#' @param max_RI_diff maximum absolute RI difference between MS/MS peaks of sample and library
#' 
#' @return data.frame with matches between features and library metabolites

#' @export
#' @importFrom stats setNames

metabolite_annotation = function(mRList = NULL, library_list = NULL, feature_spectra = NULL, rt_err_thr=1, unaccept_flag=15, accept_flag=5, suffer_flag=10, acceptable_RI = 10, acceptable_PPM_err = 10, mode=NULL, max_RI_diff=30){
  
  
  ###provide support to mRList
  feature_data <- mRList$metab_ann
  if(is.null(feature_spectra)){
    message("Checking for the presence of MS_MS_spectrum in metadata\n")
    if(any(colnames(mRList$metab_ann) == "MS_MS_spectrum")){
      message("found\n")
      features_with_MSMS <- mRList$metab_ann[!is.na(mRList$metab_ann$MS_MS_spectrum), ]
      feature_spectra <- get_spectra_list_from_vector(spectra = setNames(features_with_MSMS$MS_MS_spectrum, features_with_MSMS$Feature_ID))
      cat("MS/MS spectra available for", nrow(feature_spectra), "features.\n")
    }else{
      message("not found\n")
      stop("Impossible to proceed without feature spectra information")
    }
  }
  
  #check Retention Time similarity
  cat("Checking precursos RT...\n")
  RT = check_RT(feature_data = feature_data , reference = library_list$lib_precursor, rt_err_thr=rt_err_thr)
  
  #Calculating PPM error
  cat("Checking precursos PPM errors..\n")
  mass = check_mass(feature_data = feature_data , reference = library_list$lib_precursor, unaccept_flag=unaccept_flag, accept_flag=accept_flag, suffer_flag=suffer_flag)
  
  # Merging ideal candidate in term of Retention Time and PPM error
  cat("Finding candidate with appropriate RT and PPM error...\n")
  RT_mass = check_RT_mass (RT , mass, reference= library_list$lib_precursor )
  
  cat("Calculating RIs...\n")
  #establishing new sample data by calculating relative intensity
  RI_sample = RI_sample_data (feature_spectra= feature_spectra, accept_RI = acceptable_RI )
  
  #establishing new library data by calculating relative intensity
  #RI_lib = RI_lib_data ( reference = reference_spectra , RT_mass, acceptable_RI = acceptable_RI)
  
  #calculating and selecting proper PPM error for the most intensive peaks for mass-mass
  cat("Checking MS/MS peaks...\n")
  #intense_peak <- check_intense_peak(RT_mass = RT_mass, RI_lib = library_list$lib_peaks, reference = library_list$lib_precursor, RI_sample = RI_sample, n_peaks=n_peaks, acceptable_PPM_err = acceptable_PPM_err, mode = mode)
  
  intense_peak <- peak_matching(RT_mass = RT_mass, RI_lib = library_list$lib_peaks, reference = library_list$lib_precursor, lib_peaks_data=library_list$lib_peaks_data, mode = mode, RI_sample = RI_sample, ppm_err = acceptable_PPM_err, intensity = max_RI_diff)
  
  #intense_peak_df <- do.call(rbind, intense_peak$matched_peaks)
  
  cat("# annotated features:", length(unique(intense_peak$matched_peaks$Feature_ID)), "\n")
  cat("# metabolites:", length(unique(intense_peak$matched_peaks$ID)), "\n")
  
  return(intense_peak)
}
