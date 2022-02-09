#' Metabolite Annotation
#'
#'
#' @export
#'
metabolite_annotation = function(feature_data = NULL , reference = NULL , feature_spectra = NULL, reference_spectra= NULL,
                                 rt_window=rt_window, rt_err_thr=rt_err_thr,
                                 unaccept_flag=unaccept_flag, accept_flag=accept_flag, suffer_flag=suffer_flag,
                                 acceptable_RI = acceptable_RI,
                                 n_peaks=n_peaks, acceptable_PPM_err = acceptable_PPM_err ){


  #check Retention Time similarity
  RT = check_RT(feature_data = feature_data , reference = reference, rt_window=rt_window, rt_err_thr=rt_err_thr)

  #Calculating PPM error
  mass = check_mass(feature_data = feature_data , reference = reference, unaccept_flag=unaccept_flag, accept_flag=accept_flag, suffer_flag=suffer_flag)

  # Merging ideal candidate in term of Retention Time and PPM error
  RT_mass = check_RT_mass (RT , mass, reference= reference )

  #establishing new sample data by calculating relative intensity
  RI_sample = RI_sample_data (feature_spectra= feature_spectra, acceptable_RI = acceptable_RI )

  #establishing new library data by calculating relative intensity
  RI_lib = RI_lib_data ( reference = reference_spectra , RT_mass, acceptable_RI = acceptable_RI)

  #calculating and selecting proper PPM error for the most intensive peaks for mass-mass
  intense_peak = check_intense_peak(RT_mass, RI_lib , RI_sample, n_peaks=n_peaks, acceptable_PPM_err = acceptable_PPM_err)

  return(intense_peak)
}
