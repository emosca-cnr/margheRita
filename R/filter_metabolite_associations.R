#' Filter metabolite-feature pairs
#' @param  mRList mRList object
#' @export
#' 

filter_metabolite_associations <- function(mRList=NULL){
  
  #Level, then Mass_Status, then Peaks_found_ppm_RI, and finally RT_error.
  
  out_levels <- mRList$metabolite_identification$associations
  
  cat("Filtering...\n")
  cat("before filtering:", nrow(out_levels), "\n")
  out_levels <- unique(out_levels[, -which(colnames(out_levels) == "ID_peaks")]) #remove the redundant ID_peak
  
  ####### by LEVEL 	######
  for(lev in c(1, 2)){
    lev_m <- out_levels$ID[out_levels$Level==lev]
    if(length(lev_m) > 0){
      idx_rm <- which(out_levels$ID %in% lev_m & out_levels$Level > lev)
      if(length(idx_rm) > 0){
        out_levels <- out_levels[-idx_rm, ]
      }
    }
  }
  
  #level 3 mz,rt annotated also as level 3 mz
  lev_m <- out_levels$ID[out_levels$Level==3 & out_levels$Level_note ==  "mz, rt"]
  if(length(lev_m) > 0){
    idx_rm <- which(out_levels$ID %in% lev_m & out_levels$Level_note == "mz")
    if(length(idx_rm) > 0){
      out_levels <- out_levels[-idx_rm, ]
    }
  }
  cat("filtered:", nrow(out_levels), "\n")
  
  
  ####### by mass status 	######
  if(any(duplicated(out_levels$ID))){
    mass_status_class <- as.numeric(factor(out_levels$mass_status, levels = c("super", "acceptable", "suffer")))
    ms_num_set <- sort(unique(mass_status_class))
    if(length(ms_num_set)>1){
      for(ms in ms_num_set[-length(ms_num_set)]){ #the highest value does have higher values
        feat_selected <- out_levels$ID[mass_status_class==ms] #IDs that have level ms
        if(length(feat_selected) > 0){
          idx_rm <- which(out_levels$ID %in% feat_selected & mass_status_class > ms) #remove mass_status higher than ms for the same IDs
          if(length(idx_rm) > 0){
            out_levels <- out_levels[-idx_rm, ]
            mass_status_class <- mass_status_class[-idx_rm]
          }
        }
      }
    }
  }
  ####### by Peaks_found_ppm_RI ####
  if(any(duplicated(out_levels$ID))){
    num_peaks <- out_levels$peaks_found_ppm_RI
    num_peaks[is.na(num_peaks)] <- 0
    num_peaks_set <- unique(sort(num_peaks, decreasing = T))
    if(any(num_peaks>0)){
      for(ps in num_peaks_set[-length(num_peaks_set)]){ #the highest value does have higher values
        feat_selected <- out_levels$ID[num_peaks==ps] #IDs that have level ms
        if(length(feat_selected) > 0){
          idx_rm <- which(out_levels$ID %in% feat_selected & num_peaks > ps) #remove mass_status higher than ms for the same IDs
          if(length(idx_rm) > 0){
            out_levels <- out_levels[-idx_rm, ]
            num_peaks <- num_peaks[-idx_rm]
          }
        }
      }
    }
  }
  
 # if(any(duplicated(out_levels$ID))){
 #    ##### RT_error #####
 #    num_peaks <- out_levels$peaks_found_ppm_RI
 #    num_peaks[is.na(num_peaks)] <- 0
 #    rt_err <- out_levels$RT_err
 #    rt_err[is.na(rt_err)] <- 1000
 #    out_levels <- out_levels[order(out_levels$ID, out_levels$ppm_error, -num_peaks, rt_err), ]
 #    out_levels <- out_levels[!duplicated(out_levels$ID), ]
 #  }
  
  mRList$metabolite_identification$associations <- out_levels
  
  return(mRList)
  
}