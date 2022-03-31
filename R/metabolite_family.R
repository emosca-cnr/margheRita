#' Metabolite Family
#'
#' @export
#'
metabolite_family = function(feature_data = NULL, reference = NULL, rt_err_thr = 2, unaccept_flag= 1000, RT_mass=NULL, RI_lib=NULL , RI_sample=NULL,  lib_peaks_data=NULL, mode=c("POS", "NEG")) {
  
  
  ###RT
  RT <- vector("list", nrow(reference))
  names(RT) <- reference$ID
  
  for (k in 1:dim(reference)[1]){
    
    RT[[k]] <- data.frame(Feature_ID= feature_data$Feature_ID, RT_err=abs(reference$rt[k]- feature_data$rt), stringsAsFactors = F)
    RT[[k]]$RT_flag <- RT[[k]]$RT_err < rt_err_thr
    RT[[k]] <- RT[[k]][RT[[k]]$RT_flag, ]
  }
  
  
  ####mass
  mass = vector("list", nrow(reference))
  names(mass) = reference$ID
  
  
  for (j in 1:dim(reference)[1]) {
    
    ppm_error = abs((reference$mz[j] - feature_data$mz)) / reference$mz[j] * 1000000
    mass[[j]] = data.frame( Feature_ID = feature_data$Feature_ID , ppm_error= ppm_error,  stringsAsFactors = F)
    
    if(nrow(mass[[j]]) > 0){   #issue with empty data.frame
      
      mass[[j]]$mass_flag = mass[[j]]$ppm_error < unaccept_flag
      mass[[j]] = mass[[j]][mass[[j]]$mass_flag,]
    }
  }
  
  ###merge RT & mass
  RT_mass = check_RT_mass (RT , mass, reference= lib_data$lib_precursor)
  
  
  
  #### number of the match peaks
  
  #ensure the same order
  lib_peaks_data <- lib_peaks_data[match(names(RI_lib), lib_peaks_data$ID), ]
  
  for(z in 1:length(RT_mass)){
    RT_mass[[z]]$peaks_found_ppm_RI <- 0
    RT_mass[[z]]$peaks_found_flag <- "red"
    RT_mass[[z]]$ID <- names(RT_mass)[z]
    
    CAS_z <- reference$CAS[reference$ID == names(RT_mass)[z]]
    z_peaks_all <- RI_lib[lib_peaks_data$CAS == CAS_z & lib_peaks_data$Collision_energy == mode ]
    #z_peaks_all = z_peaks_all[sapply(z_peaks_all, function(x) length (x)[1]) > 0]
    
    
    
    if(length(z_peaks_all)>0){
      
      for(zzzz in 1:length(z_peaks_all)){
        z_peaks <- as.data.frame(z_peaks_all[[zzzz]])
        
        RT_mass[[z]]$ID_peaks <- names(z_peaks_all)[zzzz]
        
        for(zi in 1:nrow(RT_mass[[z]])){
          zi_peaks <- as.data.frame(RI_sample[names(RI_sample) == RT_mass[[z]]$Feature_ID[zi]])
          
          cat("processing zi", zi, "\n")
          
          pmm_error_matrix <- matrix(0, nrow(z_peaks), nrow(zi_peaks))
          pmm_error_matrix_flags <- pmm_error_matrix
          
          RI = pmm_error_matrix
          RI_flags = pmm_error_matrix
          
          for(i in 1:nrow(pmm_error_matrix)){
            for(j in 1:ncol(pmm_error_matrix)){
              pmm_error_matrix[i, j] = abs(z_peaks[i, 1] - zi_peaks[j, 1]) / z_peaks[i, 1] * 100
              
              RI[i, j] = abs(z_peaks[i, 2] - zi_peaks[j, 2])
            }
          }
          
          #count the peaks that are found
          pmm_error_matrix_flags[pmm_error_matrix <10] <- 1
          RI_flags[RI < 30] <- 2
          
          flags_matrix <- pmm_error_matrix_flags + RI_flags
          
          if(nrow(flags_matrix) > ncol(flags_matrix)) {
            RT_mass[[z]]$peaks_found_ppm_RI[zi] <- sum(sign(colSums(flags_matrix==3)))
          }
          
          if(nrow(flags_matrix) < ncol(flags_matrix)) {
            RT_mass[[z]]$peaks_found_ppm_RI[zi] <- sum(sign(rowSums(flags_matrix==3)))
          }
          
          
          #define the flags
          RT_mass[[z]]$peaks_found_flag[RT_mass[[z]]$peaks_found_ppm_RI > 1] <- "yellow"
          RT_mass[[z]]$peaks_found_flag[RT_mass[[z]]$peaks_found_ppm_RI > 2] <- "green"
          
          
        }
      }
    }
  }
  #filtering RT_mass by deleting the candidates with bad ppm error and RI
  for (k in 1:length(RT_mass)) {
    RT_mass[[k]] <- RT_mass[[k]][RT_mass[[k]]$peaks_found_flag != "red", ] }
  
  
  #filter the RT_mass by deleting the empty data.frame;
  RT_mass <- RT_mass[sapply(RT_mass, function(x) dim(x)[1]) > 0]
  
  RT_mass <- do.call(rbind, RT_mass)
  
  metabolites <- merge(RT_mass,lib_peaks_data, by.x="ID_peaks", by.y="ID", all.x=T, suffi)
  metabolites <- merge(metabolites, reference, by=c("ID", "CAS"))
  
  return(metabolites)
  
}
