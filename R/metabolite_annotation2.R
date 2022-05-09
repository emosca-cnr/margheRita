#' Metabolite Annotation
#' @export
#'
#'
metabolite_annotation2 <- function(reference=NULL, feature_data=NULL, RI_lib=NULL, lib_peaks_data=NULL, RI_sample=NULL, rt_err_thr=1, unaccept_flag=15, accept_flag=5, suffer_flag=10 ,mode=c("NEG","POS"),n_peaks=1, acceptable_PPM_err = 10){

  #check Retention Time similarity
  RT <- check_RT(reference=reference, feature_data=feature_data, rt_err_thr=rt_err_thr)

  #Calculating PPM error
  mass <- check_mass(reference=reference ,feature_data=feature_data, unaccept_flag=unaccept_flag, accept_flag=accept_flag, suffer_flag=suffer_flag)

  # Merging ideal candidate in term of Retention Time and PPM error
  RT_mass <- check_RT_mass (RT=RT , mass=mass, reference=reference)

  #MS/MS spectra similarity by calculating PPM error for intensive peaks
  intense_peak <- check_intense_peak(RT_mass=RT_mass, RI_lib=RI_lib, reference=reference,
                                     lib_peaks_data=lib_peaks_data, mode=mode,
                                     RI_sample=RI_sample, n_peaks=n_peaks,
                                     acceptable_PPM_err=acceptable_PPM_err)

  return(intense_peak)
}
