#' check the ppm error for the intense peaks
#'
#'
#' @export
#'
#' @param RT_mass RT_mass
#' @param RI_lib RI_lib
#' @param RI_sample RI_sample
#' @param reference reference
#' @param lib_peaks_data lib_peaks_data
#' @param n_peaks n_peaks
#' @param acceptable_PPM_err acceptable_PPM_err
#' @param mode mode
#'
check_intense_peak = function(RT_mass=NULL, RI_lib=NULL , RI_sample=NULL, reference=NULL, lib_peaks_data=NULL, n_peaks=1, acceptable_PPM_err = 10, mode=c("POS", "NEG")) {

  ans <- vector("list", length(RT_mass)) #a novel data to not touch the input
  names(ans) <- names(RT_mass)

  #ensure the same order
  lib_peaks_data <- lib_peaks_data[match(names(RI_lib), lib_peaks_data$ID), ]

  #z cycles trhough the precursor library entries
  for(z in 1:length(RT_mass)){
    #cat(z)
    #ans[[z]]$intense_peak <- F
    #ans[[z]]$intense_peak2 <- NA
    #ans[[z]]$intense_peak3 <- NA

    #clean the library to keep only metabolites that appear in the sample
    #z_peaks <- as.data.frame(RI_lib[names(RI_lib) == names(RT_mass)[z]])
    CAS_z <- reference$CAS[reference$ID == names(RT_mass)[z]]
    z_peaks_all <- RI_lib[lib_peaks_data$CAS == CAS_z & lib_peaks_data$Collision_energy == mode ]

    if(length(z_peaks_all)>0){
      #cycle thourgh the peaks that are available for the given CAS at precursor level
      for(zzzz in 1:length(z_peaks_all)){

        z_peaks <- z_peaks_all[[zzzz]]

        #cycle thourgh the candidates
        for(zi in 1:nrow(RT_mass[[z]])){

          zi_peaks <- as.data.frame(RI_sample[names(RI_sample) == RT_mass[[z]]$Feature_ID[zi]])


          #the first
          table_values_a <- sort(unique(z_peaks[, 2]), decreasing = T) #library peaks
          table_values_b <- sort(unique(zi_peaks[, 2]), decreasing = T) #sample peaks

          temp_flag_peaks <- rep(FALSE, n_peaks)
          for(ii in 1:min(n_peaks, length(table_values_a))){

            a = z_peaks[z_peaks[,2] == table_values_a[ii], , drop=FALSE]
            if(length(table_values_b) < ii){
              break
            }
            b = zi_peaks[zi_peaks[,2] == table_values_b[ii], , drop=FALSE]

            PPM_err <- calc_ppm_err(a, b)

            #RT_mass[[z]][zi, 6+ii] <- any(PPM_err < acceptable_PPM_err)
            temp_flag_peaks[ii] <- any(PPM_err < acceptable_PPM_err) #

            #if(!RT_mass[[z]][zi, 6+ii]){
            #  break
            #}

            if(!temp_flag_peaks[ii]){ #if FALSE avoid the comparison for the next peaks
              break
            }

          }

          #add information about the peak to the precursor
          if(is.null(ans[[z]])){
            ans[[z]] <- data.frame(Feature_ID=RT_mass[[z]]$Feature_ID[zi], ID_peaks=names(z_peaks_all)[zzzz], intense_peaks=all(temp_flag_peaks), stringsAsFactors = F)
          }else{
            ans[[z]] <- rbind(ans[[z]], data.frame(Feature_ID=RT_mass[[z]]$Feature_ID[zi], ID_peaks=names(z_peaks_all)[zzzz], intense_peaks=all(temp_flag_peaks), stringsAsFactors = F))

          }

        } #end cycle through candidates

      } #end cycle through the peaks for that CAS

      #RT_mass[[z]] = RT_mass[[z]][apply(RT_mass[[z]][, c("intense_peak", "intense_peak2", "intense_peak3")], 1, all, na.rm=TRUE), ]
    }

    if(!is.null(ans[[z]])){
      #keep only peaks with TRUE flags
      ans[[z]] <- merge(RT_mass[[z]], ans[[z]], by="Feature_ID", sort=F)
      ans[[z]] <- data.frame(ID=names(ans)[z], ans[[z]], stringsAsFactors=F)
    }

  }

  #filter the RT_mass by deleting the empty data.frame;
  ans <- ans[sapply(ans, function(x) !is.null(x))]

  if(length(ans)>0){
    ans <- lapply(ans, function(x) x[x$intense_peaks, ])
  }else{
    message("none of the candidates matches the library peaks with current parameters\n")
    return()
  }

  if(length(ans)>0){
    ans <- ans[sapply(ans, function(x) nrow(x)>0)]

  }else{
    message("none of the candidates matches the library peaks with current parameters\n")
    return()
  }

  if(length(ans)>0){
    ans <- do.call(rbind, ans)
    ans <- merge(ans, lib_peaks_data, by.x="ID_peaks", by.y="ID", all.x=T)
    colnames(ans)[colnames(ans) == "Name"] <- "Name_peaks"
    ans <- merge(reference, ans, by=c("ID", "CAS"), all.y=T)
    ans <- ans[, c("Feature_ID", "ID", "CAS", "PubChemCID", "Name", "rt", "RT_err", "RT_flag", "mz", "ppm_error", "mass_status", "mass_flag", "ID_peaks", "Name_peaks", "intense_peaks", "Collision_energy")]
  }else{
    message("none of the candidates matches the library peaks with current parameters\n")
    return()
  }

  return(ans)

}

