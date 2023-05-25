#' Metabolite Identification
#' Run the metabolite identification pipeline
#' @param mRList mRList object
#' @param library_list library data obtained from margherita_library() function
#' @param features the feaatures to be annotated. If NULL all features in the mRList will be annotated
#' @param rt_err threshold over the retention time error
#' @param rt_best_thr Features with RT error < rt_best_thr  will be labelled as "super", those with rt_best_thr < RT error <= rt_err_thr as "acceptable", the remaining as "unacceptable".
#' @param accept_flag A number with default value of 5. PPM errors < accept_flag will be tagged as "super", while those > accept_flag and < suffer_flag will be tagged as "acceptable"
#' @param suffer_flag A number with default value of 10. PM errors above this value and < unaccept_flag will be tagged as "suffer"
#' @param unaccept_flag A number with default value of 15. The maximum PPM error must be less than this value. and those above this number will be eliminated.
#' @param min_RI numeric parameter. the default value is 10. it is a maximum relative intensity that is kept in sample. since low intense 
#' peaks could be noise, it is filtering sample dataset by deleting the relative intensity lower then accept_RI.
#' @param ppm_err A number with default value of 10. The maximum PPM error must be less than this value. and those above this number will be eliminated.
#' @param RI_err maximum absolute RI difference between MS/MS peaks of sample and library
#' @param RI_err_type type of RI error calculation.
#' @param filter whether to filter metabolite-feature associations or not.
#' @param lib_ann_fields columns of library_list$lib_precursor that will be added to metabolite annotations
#' @return data.frame with matches between features and library metabolites
#' @export
#' @importFrom stats setNames

metabolite_identification <- function(mRList = NULL, features = NULL, library_list = NULL, rt_err=1, rt_best_thr=0.5, unaccept_flag=20, accept_flag=5, suffer_flag=10, min_RI = 10, ppm_err = 20, RI_err=20, RI_err_type="rel", filter=FALSE, lib_ann_fields=c("ID", "Name", "PubChemCID")){
  
  if(is.null(features)){
    idx <- 1:nrow(mRList$metab_ann) 
  }else{
    idx <- which(mRList$metab_ann$Feature_ID %in% features) 
  }
  
  cat("# Features considered: ", length(idx), ".\n")
  
  #extract feature information from metabolite metadata
  out_levels <- unique(mRList$metab_ann[idx, c("Feature_ID", "rt", "mz", "MS_MS_spectrum")])
  
  
  ### check retention time and mass
  if(length(library_list$lib_precursor$rt)>0){
    cat("Checking RT... ")
    RT <-  check_RT(feature_data = out_levels, reference = library_list$lib_precursor, rt_err_thr=rt_err, rt_best_thr = rt_best_thr)
    cat("done\n")
  }
  
  cat("Checking precursor mz... ")
  mass <- check_mass(feature_data = out_levels, reference = library_list$lib_precursor, unaccept_flag=unaccept_flag, accept_flag=accept_flag, suffer_flag=suffer_flag)
  cat("done\n")
  
  if(length(library_list$lib_precursor$rt)>0){
    cat("Integrating RT and mz results... ")
    RT_mass <- check_RT_mass(RT = RT, mass = mass)
    cat("done\n")
  }
  
  cat("Assembling output for precursors...")
  ### features with appropriate mass
  mass_ok <- stack(lapply(mass, function(x) x$Feature_ID))
  mass_ok <- merge(library_list$lib_precursor, mass_ok, by.x="ID", by.y="ind", all.y=T)
  colnames(mass_ok) <- replace(colnames(mass_ok), list = colnames(mass_ok)=="values", "Feauture_ID")
  
  out_levels <- merge(out_levels, data.frame(ID=rep(names(mass), unlist(lapply(mass, nrow))), do.call(rbind, lapply(mass, function(x) x[, colnames(x) != "mz"])), stringsAsFactors = F), by="Feature_ID", all=T)
  out_levels <- merge(library_list$lib_precursor, out_levels, by="ID", all = T, suffixes = c("_lib", ""))
  
  ### features with appropriate mass and RT
  if(length(library_list$lib_precursor$rt)>0){
    
    RT_mass_ok <- stack(lapply(RT_mass, function(x) x$Feature_ID))
    RT_mass_ok <- merge(library_list$lib_precursor, RT_mass_ok, by.x="ID", by.y="ind", all.y=T)
    
    out_levels <- merge(out_levels, data.frame(ID=rep(names(RT_mass), unlist(lapply(RT_mass, nrow))), do.call(rbind, lapply(RT_mass, function(x) x[, c("Feature_ID", "rt", "RT_err", "RT_class", "RT_flag")])), stringsAsFactors = F), by=c("ID", "Feature_ID", "rt"), all=T)
    
  }else{
    out_levels$RT_err <- NA
    out_levels$RT_class <- NA
    out_levels$RT_flag <- NA
    out_levels$rt_lib <- NA
  }
  
  cat("done\n")
  
  
  ###MS/MS matching - features that match at mass level and have MSMSspectra
  cat("MS/MS analysis...\n")
  feature_spectra_list <- unique(out_levels$Feature_ID[which(!is.na(out_levels$MS_MS_spectrum) & out_levels$mass_flag)]) ## features already associated at least with mass and MS/MS available
  
  feature_spectra_list <- get_spectra_list_from_vector(spectra = setNames(mRList$metab_ann$MS_MS_spectrum[mRList$metab_ann$Feature_ID %in% feature_spectra_list], mRList$metab_ann$Feature_ID[mRList$metab_ann$Feature_ID %in% feature_spectra_list])) #list of features and their spectra
  
  RI_sample <- RI_sample_data(feature_spectra=feature_spectra_list, accept_RI = min_RI) #define RI
  
  #remove from mass those features that do not appear in RI_sample, because do not have MSMS spectra
  mass <- lapply(mass, function(x) x[x$Feature_ID %in% names(RI_sample), ])
  mass <- mass[unlist(lapply(mass, nrow)) > 0]
  
  intense_peak <- peak_matching(library_list = library_list, library_matched_features = mass, RI_sample = RI_sample, ppm_err = ppm_err, intensity = RI_err, RI_diff_type = RI_err_type)
  cat("done\n")
  
  
  cat("Assembling final output... ")
  
  intense_peak_unq <- unique(intense_peak$matched_peaks[, c("ID", "Feature_ID")])
  intense_peak_unq <- merge(library_list$lib_precursor, intense_peak_unq, by="ID", all.y=T)
  
  
  out_levels <- merge(out_levels, intense_peak$matched_peaks[, c("ID", "Feature_ID", "ID_peaks", "peaks_found_ppm_RI", "matched_peaks_ratio", "precursor_in_MSMS")], by=c("ID", "Feature_ID"), all=T)
  
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
  out_levels$Level_note[idx] <- ""
  
  idx <- which(!is.na(out_levels$peaks_found_ppm_RI) & out_levels$mass_flag & out_levels$RT_flag)
  out_levels$Level[idx] <- 1
  out_levels$Level_note[idx] <- ""
  
  
  ### overall Reliability
  #out_levels$Reliability <- ""
  #out_levels$Reliability
  
  ### summary
  out_levels <- out_levels[out_levels$Level != "", ]
  
  #out_levels <- out_levels[, c("Feature_ID", "rt", "mz", "MS_MS_spectrum", "RT_err", "RT_flag", "ppm_error", "mass_flag", "mass_status", "ID_peaks", "peaks_found_ppm_RI", "precursor_in_MSMS", "ID", "Name", "rt_lib", "mz_lib", "Level", "Level_note", "CAS", "PubChemCID")]
  out_levels <- out_levels[, c("Feature_ID", "rt", "mz", "RT_err", "RT_class", "RT_flag", "ppm_error", "mass_flag", "mass_status", "ID_peaks", "peaks_found_ppm_RI", "matched_peaks_ratio", "precursor_in_MSMS", "ID", "Name", "rt_lib", "mz_lib", "Level", "Level_note")]
  
  #add the results to mRList
  mRList$metabolite_identification <- list(associations=out_levels, RI_sample=RI_sample, MS_MS_info=intense_peak$matrices)
  
  ## FILTERING
  #metabolite-feature association of level 1 imply that these do not appear in other pairs, and so on for other levels
  if(filter){
    mRList <- filter_metabolite_associations(mRList)
  }
  
  
  ### summary associations and Feature ID assignment
  #summary associations
  mRList$metabolite_identification$associations_summary <- unique(out_levels[, c("Feature_ID", "ID", "Name", "Level", "Level_note")])
  
  cat("Adding annotations to metab_ann...\n")
  #add name to metab_ann
  feat_name <- split(mRList$metabolite_identification$associations_summary$ID, mRList$metabolite_identification$associations_summary$Feature_ID)
  
  for(i_col in lib_ann_fields){
    mapped_feat <- as.data.frame(unlist(lapply(feat_name, function(x) paste0(sort(library_list$lib_precursor[match(x, library_list$lib_precursor$ID), i_col]), collapse = ";"))))
    colnames(mapped_feat) <- i_col
    mRList$metab_ann <- merge(mRList$metab_ann, mapped_feat, by.x="Feature_ID", by.y=0, all.x = T, sort = F)
  }

  mRList$metab_ann <- mRList$metab_ann[order(mRList$metab_ann$Feature_ID), ]
  
  cat("done\n")
  
  return(mRList)
  
}
