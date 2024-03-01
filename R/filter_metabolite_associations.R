#' Filter metabolite-feature pairs
#' @param  mRList mRList object
#' @export
#' 

filter_metabolite_associations <- function(mRList=NULL){
  
  out_levels <- mRList$metabolite_identification$associations

  cat("before filtering:", nrow(out_levels), "\n")
  
  cat("Filtering metabolites-features 1-to-many associations...\n")
  
  cat("Filtering associations with MS/MS... ")

  out_levels$ID_Feature_ID <- apply(out_levels[, c("ID", "Feature_ID")], 1, paste0, collapse="_")
  out_levels$ID_Feature_ID <- gsub(" +", "", out_levels$ID_Feature_ID)
  
  out_levels_NA <- out_levels[is.na(out_levels$peaks_found_ppm_RI), ]
  out_levels_tofilter <- out_levels[!is.na(out_levels$peaks_found_ppm_RI), ]
  
  ### To resolve the one-to-many associations between (f_i,m_j) pairs and MS/MS spectra s_jk, for each pair (f_i, m_j), select the s_jk that has:
  ### - highest number of peaks;
  ### - highest peaks ratio;

  out_levels_tofilter <- out_levels_tofilter[order(-out_levels_tofilter$peaks_found_ppm_RI, -out_levels_tofilter$matched_peaks_ratio), ]
  
  idx_rm <- which(duplicated(out_levels_tofilter$ID_Feature_ID))
  if(length(idx_rm) > 0){
    #cat("Removing", length(idx_rm), "\n")
    out_levels_tofilter <- out_levels_tofilter[-idx_rm, ]
  }
  
  #To resolve the one-to-many associations between m_j and f_i, for each m_j, select the f_i that has:
  out_levels_tofilter <- out_levels_tofilter[order(-out_levels_tofilter$matched_peaks_ratio, out_levels_tofilter$Level, out_levels_tofilter$ppm_error), ]
  
  idx_rm <- which(duplicated(out_levels_tofilter$ID))
  if(length(idx_rm) > 0){
    #cat("Removing", length(idx_rm), "\n")
    out_levels_tofilter <- out_levels_tofilter[-idx_rm, ]
  }
  out_levels <- rbind(out_levels_tofilter, out_levels_NA)
  out_levels$ID_Feature_ID <- NULL
  
  cat(nrow(out_levels), "\n")
  
  #Level, then Mass_Status, then Peaks_found_ppm_RI, and finally RT_error.
  
  cat("Filtering associations without MS/MS...\n")
  
  ####### by LEVEL 	######
  cat("\tby level... ")
  for(lev in c(1, 2)){
    lev_m <- out_levels$ID[out_levels$Level==lev]
    if(length(lev_m) > 0){
      idx_rm <- which(out_levels$ID %in% lev_m & out_levels$Level > lev)
      if(length(idx_rm) > 0){
        #cat(unique(out_levels$Level[idx_rm]), "\n")
        out_levels <- out_levels[-idx_rm, ]
      }
    }
  }
  
  lev_m <- out_levels$ID[out_levels$Level==3 & out_levels$Level_note ==  "mz, rt"]
  if(length(lev_m) > 0){
    idx_rm <- which(out_levels$ID %in% lev_m & out_levels$Level_note == "mz")
    if(length(idx_rm) > 0){
      #cat(unique(out_levels$Level[idx_rm]), "\n")
      out_levels <- out_levels[-idx_rm, ]
    }
  }
  cat(nrow(out_levels), "\n")
  
  ####### by mass status 	######
  if(any(duplicated(out_levels$ID))){
    cat("\tby mass... ")
    mass_status_class <- as.numeric(factor(out_levels$mass_status, levels = c("super", "acceptable", "suffer")))
    ms_num_set <- sort(unique(mass_status_class))
    if(length(ms_num_set)>1){
      for(ms in ms_num_set[-length(ms_num_set)]){ #the highest value does have higher values
        feat_selected <- out_levels$ID[mass_status_class==ms] #IDs that have level ms
        if(length(feat_selected) > 0){
          idx_rm <- which(out_levels$ID %in% feat_selected & mass_status_class > ms) #remove mass_status higher than ms for the same IDs
          if(length(idx_rm) > 0){
            #cat(unique(out_levels$Level[idx_rm]), "\n")
            out_levels <- out_levels[-idx_rm, ]
            mass_status_class <- mass_status_class[-idx_rm]
          }
        }
      }
    }
  }
  cat(nrow(out_levels), "\n")
  
  ##### RT_class #####
  if(any(duplicated(out_levels$ID))){
    
    cat("\tby RT... ")
    
    rt_class <- as.numeric(factor(out_levels$RT_class, levels = c("super", "acceptable")))
    rt_class_set <- sort(unique(rt_class))
    
    if(length(rt_class_set)>1){
      for(rt_i in rt_class_set[-length(rt_class_set)]){ #the highest value does have higher values
        feat_selected <- out_levels$ID[rt_class==rt_i] #IDs that have level ms
        if(length(feat_selected) > 0){
          idx_rm <- which(out_levels$ID %in% feat_selected & rt_class > rt_i) #remove mass_status higher than ms for the same IDs
          if(length(idx_rm) > 0){
            #cat(unique(out_levels$Level[idx_rm]), "\n")
            out_levels <- out_levels[-idx_rm, ]
            rt_class <- rt_class[-idx_rm]
          }
        }
      }
    }
  }
  cat(nrow(out_levels), "\n")
  
  
  #### filter over features, to remove the worst metabolite IDs
  cat("Filtering features-metabolites 1-to-many associations...\n")
  
  ####### by LEVEL 	######
  cat("\tby level... ")
  for(lev in c(1, 2)){
    lev_m <- out_levels$Feature_ID[out_levels$Level==lev]
    if(length(lev_m) > 0){
      idx_rm <- which(out_levels$Feature_ID %in% lev_m & out_levels$Level > lev)
      if(length(idx_rm) > 0){
        #cat(unique(out_levels$Level[idx_rm]), "\n")
        out_levels <- out_levels[-idx_rm, ]
      }
    }
  }
  
  #level 3 mz,rt annotated also as level 3 mz
  lev_m <- out_levels$Feature_ID[out_levels$Level==3 & out_levels$Level_note ==  "mz, rt"]
  if(length(lev_m) > 0){
    idx_rm <- which(out_levels$Feature_ID %in% lev_m & out_levels$Level_note == "mz")
    if(length(idx_rm) > 0){
      #cat(unique(out_levels$Level[idx_rm]), "\n")
      out_levels <- out_levels[-idx_rm, ]
    }
  }
  cat(nrow(out_levels), "\n")
  
  
  ####### by mass status 	######
  if(any(duplicated(out_levels$Feature_ID))){
    cat("\tby mass... ")
    mass_status_class <- as.numeric(factor(out_levels$mass_status, levels = c("super", "acceptable", "suffer")))
    ms_num_set <- sort(unique(mass_status_class))
    if(length(ms_num_set)>1){
      for(ms in ms_num_set[-length(ms_num_set)]){ #the highest value does have higher values
        feat_selected <- out_levels$Feature_ID[mass_status_class==ms] #IDs that have level ms
        if(length(feat_selected) > 0){
          idx_rm <- which(out_levels$Feature_ID %in% feat_selected & mass_status_class > ms) #remove mass_status higher than ms for the same IDs
          if(length(idx_rm) > 0){
            #cat(unique(out_levels$Level[idx_rm]), "\n")
            out_levels <- out_levels[-idx_rm, ]
            mass_status_class <- mass_status_class[-idx_rm]
          }
        }
      }
    }
  }
  cat(nrow(out_levels), "\n")
  
  ####### by Peaks_found_ppm_RI ####
  if(any(duplicated(out_levels$Feature_ID))){
    cat("\tby number of matching MS/MS peaks... ")
    num_peaks <- out_levels$peaks_found_ppm_RI
    num_peaks[is.na(num_peaks)] <- 0
    num_peaks_set <- unique(sort(num_peaks, decreasing = T))
    if(any(num_peaks>0)){
      for(ps in num_peaks_set[-length(num_peaks_set)]){ #the highest value does have less peaks
        feat_selected <- out_levels$Feature_ID[num_peaks==ps] #IDs with ps peaks
        if(length(feat_selected) > 0){
          idx_rm <- which(out_levels$Feature_ID %in% feat_selected & num_peaks < ps) #remove mass_status higher than ms for the same IDs >
          if(length(idx_rm) > 0){
            #cat(unique(out_levels$Level[idx_rm]), "\n")
            out_levels <- out_levels[-idx_rm, ]
            num_peaks <- num_peaks[-idx_rm]
          }
        }
      }
    }
  }
  cat(nrow(out_levels), "\n")
  
  if(any(duplicated(out_levels$Feature_ID))){
    
    ##### RT_class #####
    cat("\tby RT... ")
    rt_class <- as.numeric(factor(out_levels$RT_class, levels = c("super", "acceptable")))
    rt_class_set <- sort(unique(rt_class))
    
    if(length(rt_class_set)>1){
      for(rt_i in rt_class_set[-length(rt_class_set)]){ #the highest value does have higher values
        feat_selected <- out_levels$Feature_ID[rt_class==rt_i] #IDs that have level ms
        if(length(feat_selected) > 0){
          idx_rm <- which(out_levels$Feature_ID %in% feat_selected & rt_class > rt_i) #remove mass_status higher than ms for the same IDs
          if(length(idx_rm) > 0){
            #cat(unique(out_levels$Level[idx_rm]), "\n")
            out_levels <- out_levels[-idx_rm, ]
            rt_class <- rt_class[-idx_rm]
          }
        }
      }
    }
  }
  cat(nrow(out_levels), "\n")
  
  mRList$metabolite_identification$associations <- out_levels
  
  mRList$metabolite_identification$associations_summary <- unique(out_levels[, c("Feature_ID", "ID", "Name", "Level", "Level_note")])
  
  return(mRList)
  
}