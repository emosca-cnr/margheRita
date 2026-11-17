#' MS/MS Spectra Visualization
#' Visualize the spectrum of a metabolite and each of the features associted with it.
#'
#' @param mRList mRList object that contains the results of metabolite identification
#' @param mR_library library list used during annotation
#' @param metabolite_id ID of the metabolite to be visualized
#' @param peak_number if TRUE write the peak number
#' @param type mirrored or overlapped visualization
#' @returns list of ggplot objects
#' @export

visualize_associated_spectra <- function(mRList=NULL, mR_library=NULL, metabolite_id=NULL, type=c("mirrored", "overlapped"), peak_number=FALSE){

  type <- match.arg(type)

  metabolite_name <- unique(mRList$metabolite_identification$associations$Name[mRList$metabolite_identification$associations$ID == metabolite_id])

  if(length(metabolite_name) == 0){
    stop("Can't find ", metabolite_id, " between level 1-2 associations.\n")
  }
  if(!any(mRList$metabolite_identification$associations$Level[mRList$metabolite_identification$associations$Name %in% metabolite_name] %in% c(1, 2))){
    stop("Can't find ", metabolite_id, " between level 1-2 associations.\n")
  }

  #get the features associated with metabolite_id
  associated_features <- unique(mRList$metabolite_identification$associations[mRList$metabolite_identification$associations$ID == metabolite_id & mRList$metabolite_identification$associations$Level %in% c(1, 2), c("Feature_ID", "ID_peaks"), drop=F])

  RI_sample <- mRList$metabolite_identification$RI_sample
  RI_lib <- mR_library$lib_peaks

  p_list <- vector("list", length = nrow(associated_features))
  for(i in 1:nrow(associated_features)){

    feature_id <- associated_features$Feature_ID[i]
    peak_id <- associated_features$ID_peaks[i]

    feature_spectra <- RI_sample[[as.character(feature_id)]]
    library_spectra <- RI_lib[[as.character(peak_id)]]
    match_matrix <- mRList$metabolite_identification$MS_MS_info[[metabolite_id]][[peak_id]][[as.character(feature_id)]]$flags

    
    p_list[[i]] <- visualize_feature_metabolite_spectra(spectra_top = library_spectra, spectra_bottom = feature_spectra, name_top = peak_id, name_bottom = feature_id, type = type, peak_number = peak_number, match_matrix = match_matrix)
    

  }
  
  return(p_list)

}
