#' Filter metabolite-feature pairs
#' @param  mRList mRList object
#' @export
#' @return filtered mRList object


filter_metabolite_associations <- function(mRList=NULL){
  
  out_levels <- mRList$metabolite_identification$associations
  
  cat("Initial associations:", nrow(out_levels), "\n")
  
  out_levels$ID_Feature_ID <- apply(out_levels[, c("ID", "Feature_ID")], 1, paste0, collapse="_")
  out_levels$ID_Feature_ID <- gsub(" +", "", out_levels$ID_Feature_ID)
  
  
  ### 1) To resolve the one-to-many associations between (f_i,m_j) pairs and MS/MS spectra s_jk, for each pair (f_i, m_j), select the s_jk that has:
  ### - highest number of peaks;
  ### - highest peaks ratio;
  cat("1) Filtering redundant MS/MS spectra for the same metabolite-feature...\n")
  cat("\tby highest number of peaks\n")
  cat("\tby highest peaks ratio...")
  
  out_levels_NA <- out_levels[is.na(out_levels$peaks_found_ppm_RI), ]
  out_levels_tofilter <- out_levels[!is.na(out_levels$peaks_found_ppm_RI), ]
  
  out_levels_tofilter <- out_levels_tofilter[order(-out_levels_tofilter$peaks_found_ppm_RI, -out_levels_tofilter$matched_peaks_ratio), ]
  
  idx_rm <- which(duplicated(out_levels_tofilter$ID_Feature_ID))
  if(length(idx_rm) > 0){
    #cat("Removing", length(idx_rm), "\n")
    out_levels_tofilter <- out_levels_tofilter[-idx_rm, ]
  }
  cat(nrow(out_levels_NA) + nrow(out_levels_tofilter), "\n")
  
  #2) To resolve the one-to-many associations between m_j and f_i, for each m_j, select the f_i that has:
  #a.	the highest ratio between #{feature MS/MS peaks that match metabolite MS/MS peaks} and #{metabolite MS/MS peaks};
  #b.	the best annotation level;
  #c.	the lowest ppm error
  cat("2) 1:n associations between metabolites and features...\n")
  cat("\tby highest peaks ratio\n")
  cat("\tby best annotaiton level\n")
  cat("\tby lowest ppm error...")
  
  out_levels_tofilter <- out_levels_tofilter[order(-out_levels_tofilter$matched_peaks_ratio, out_levels_tofilter$Level, out_levels_tofilter$ppm_error), ]
  
  idx_rm <- which(duplicated(out_levels_tofilter$ID))
  if(length(idx_rm) > 0){
    #cat("Removing", length(idx_rm), "\n")
    out_levels_tofilter <- out_levels_tofilter[-idx_rm, ]
  }
  out_levels <- rbind(out_levels_tofilter, out_levels_NA)
  out_levels$ID_Feature_ID <- NULL
  
  cat(nrow(out_levels), "\n")
  
  #	3 To resolve the one-to-many associations between m_j and f_i not supported by MS/MS spectra, for each m_j, select the f_i that has: 
  # the best annotation level;
  # the best mass match;
  # the best retention time match
  # Lowest rt error
  # Lowest ppm error
  
  cat("3) 1:n associations between metabolites and features, not supported by MS/MS...\n")
  # cat("\tby best annotaiton level\n")
  # cat("\tby best ppm match\n")
  # cat("\tby best RT match\n")
  # cat("\tby lowest RT error\n")
  # cat("\tby lowest ppm error\n")
  
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
    cat("\tby mass (category)... ")
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
    
    cat("\tby RT (category)... ")
    
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
  
  #duplicated IDs
  m_id <- unique(out_levels$ID[duplicated(out_levels$ID)])
  
  out_levels$selected <- TRUE
  
  for(i in 1:length(m_id)){
    
    #by lowest RT error
    if(sum(out_levels$selected[out_levels$ID == m_id[i]])>1){
      
      rt_err_min <- min(out_levels$RT_err[out_levels$ID == m_id[i] & out_levels$selected])
      
      if(any(out_levels$RT_class[out_levels$ID ==  m_id[i] & out_levels$selected] != "unacceptable")){
        
        if(any(abs(out_levels$RT_err[out_levels$ID ==  m_id[i] & out_levels$selected] - rt_err_min) > 1e-4)){
          
          cat("\tby RT... ")
          idx_rm <- which(out_levels$ID == m_id[i] & out_levels$selected & abs(out_levels$RT_err - rt_err_min) > 1e-4)
          out_levels$selected[idx_rm] <- FALSE
          #out_levels <- out_levels[-idx_rm, ]
          #print(out_levels[out_levels$ID == m_id[i], ])
          cat(sum(out_levels$selected), "\n")
          
        }
        
      }
    }
    
    #by lowest ppm error
    if(sum(out_levels$selected[out_levels$ID == m_id[i]])>1){
      
      ppm_err_min <- min(out_levels$ppm_error[out_levels$ID == m_id[i] & out_levels$selected])
      
      if(any(abs(out_levels$ppm_error[out_levels$ID ==  m_id[i] & out_levels$selected] - ppm_err_min) > 1e-4)){
        
        cat("\tby ppm... ")
        idx_rm <- which(out_levels$ID == m_id[i] & out_levels$selected & abs(out_levels$ppm_error - ppm_err_min) > 1e-4)
        out_levels$selected[idx_rm] <- FALSE
        #out_levels <- out_levels[-idx_rm, ]
        #print(out_levels[out_levels$ID == m_id[i], ])
        cat(sum(out_levels$selected), "\n")
        
      }
    }
    
    
  }
  
  out_levels <- out_levels[out_levels$selected, ]
  out_levels$selected <- NULL
  
  #cat(nrow(out_levels), "\n")
  
  # To resolve the one-to-many associations between f_i and m_j, for each f_i keep the m_j that has:
  #   the best level;
  # the best mass match;
  # the highest number of matching MS/MS peaks;
  # the best retention time.
  # the highest ratio between #{feature MS/MS peaks that match metabolite MS/MS peaks} and #{metabolite MS/MS peaks}
  # Lowest rt error
  # Lowest ppm error
  
  cat("4) 1:n associations between features and metabolites...\n")
  # cat("\tby best annotaiton level\n")
  # cat("\tby best ppm match\n")
  # cat("\tby highest number of matching MS/MS peaks\n")
  # cat("\tby best RT match\n")
  # cat("\tby highest ratio of MS/MS peaks\n")
  # cat("\tby lowest RT error\n")
  # cat("\tby lowest ppm error\n")
  
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
    cat("\tby mass (category)... ")
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
    cat("\tby RT (category)... ")
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
  
  ### further FILTERING
  f_id <- unique(out_levels$Feature_ID[duplicated(out_levels$Feature_ID)])
  
  out_levels$selected <- TRUE
  
  for(i in 1:length(f_id)){
    
    
    if(sum(out_levels$selected[out_levels$Feature_ID == f_id[i]])>1){
      
      ### by matched_peaks_ratio
      
      if(any(!is.na(out_levels$matched_peaks_ratio[out_levels$Feature_ID == f_id[i]]))){
        pr_max <- max(out_levels$matched_peaks_ratio[out_levels$Feature_ID == f_id[i]], na.rm = T)
        
        if(!is.na(pr_max)){
          
          if(any(abs(out_levels$matched_peaks_ratio[out_levels$Feature_ID ==  f_id[i]] - pr_max) > 1e-4)){
            
            cat("\tby matched peak ratio... ")
            idx_rm <- which(out_levels$Feature_ID == f_id[i] & abs(out_levels$matched_peaks_ratio - pr_max) > 1e-4)
            out_levels$selected[idx_rm] <- FALSE
            cat(sum(out_levels$selected), "\n")
            
          }
          
        }
      }
      
    }
    
    if(sum(out_levels$selected[out_levels$Feature_ID == f_id[i]])>1){
      
      rt_err_min <- min(out_levels$RT_err[out_levels$Feature_ID == f_id[i] & out_levels$selected])
      
      if(any(out_levels$RT_class[out_levels$Feature_ID ==  f_id[i] & out_levels$selected] != "unacceptable")){
        
        if(any(abs(out_levels$RT_err[out_levels$Feature_ID ==  f_id[i] & out_levels$selected] - rt_err_min) > 1e-4)){
          
          cat("\tby RT... ")
          idx_rm <- which(out_levels$Feature_ID == f_id[i] & out_levels$selected & abs(out_levels$RT_err - rt_err_min) > 1e-4)
          out_levels$selected[idx_rm] <- FALSE
          #out_levels <- out_levels[-idx_rm, ]
          #print(out_levels[out_levels$Feature_ID == f_id[i], ])
          cat(sum(out_levels$selected), "\n")
          
        }
        
      }
    }
    
    
    if(sum(out_levels$selected[out_levels$Feature_ID == f_id[i]])>1){
      
      
      ppm_err_min <- min(out_levels$ppm_error[out_levels$Feature_ID == f_id[i] & out_levels$selected])
      
      if(any(abs(out_levels$ppm_error[out_levels$Feature_ID ==  f_id[i] & out_levels$selected] - ppm_err_min) > 1e-4)){
        
        cat("\tby ppm... ")
        idx_rm <- which(out_levels$Feature_ID == f_id[i] & out_levels$selected & abs(out_levels$ppm_error - ppm_err_min) > 1e-4)
        out_levels$selected[idx_rm] <- FALSE
        #out_levels <- out_levels[-idx_rm, ]
        #print(out_levels[out_levels$Feature_ID == f_id[i], ])
        cat(sum(out_levels$selected), "\n")
        
      }
    }
    
    
  }
  
  out_levels <- out_levels[out_levels$selected, ]
  out_levels$selected <- NULL
  
  mRList$metabolite_identification$associations <- out_levels
  
  mRList$metabolite_identification$associations_summary <- unique(out_levels[, c("Feature_ID", "ID", "Name", "Level", "Level_note")])
  mRList$metabolite_identification$associations_summary <- unique(out_levels[, c("Feature_ID", "ID", "Name", "Level", "Level_note", "RT_err", "ppm_error", "peaks_found_ppm_RI", "matched_peaks_ratio")])
  
  
  return(mRList)
  
}




