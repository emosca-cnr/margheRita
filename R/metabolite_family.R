#' family of metabolites
#' @export
#'
metabolite_family=function(reference=NULL, feature_data=NULL, RI_lib=NULL, RI_sample=NULL,lib_peaks_data=NULL, rt_err_thr=1, unaccept_flag= 10, mode=c("NEG","POS"), ppm_err=10, intensity=30){

RT <- check_RT(reference=reference, feature_data=feature_data, rt_err_thr=rt_err_thr)

mass <- check_mass2(reference=reference, feature_data=feature_data ,unaccept_flag=unaccept_flag)

RT_mass <- check_RT_mass(RT=RT , mass=mass , reference=reference)

match_peaks = peak_matching(feature_data=feature_data ,reference=reference ,RT_mass=RT_mass,
                             RI_lib=RI_lib ,RI_sample=RI_sample ,lib_peaks_data=lib_peaks_data,
                             mode=mode, ppm_err=ppm_err, intensity=intensity)

return(match_peaks)

}
