#' metabolite_annotation
#' map between metabolite name and various stantard identifiers
#' @param mRList mRList object
#' @param reference data.frame with columns "rt", "mz", "Name", "MS_MS_spectrum"
#' @export
#' @return the annotation list, which has three elements... 1) named list of data.frames; each element of the list corresponds to an element of reference and each data.frame containes the candidates
metabolite_annotation <- function(mRList = NULL , reference = NULL){

  ## feature data extraction from mRList
  feature_data <- mRList$metab_ann
  feature_data_spectra <- setNames(mRList$metab_ann$MS_MS_spectrum, mRList$metab_ann$Feature_ID)
  feature_data_spectra <- get_spectra_list_from_vector(feature_data_spectra)

  ## reference data extraction from library
  reference_spectra <- setNames(reference$MS_MS_spectrum, reference$Name)
  reference_spectra <- get_spectra_list_from_vector(reference_spectra)
  reference <- reference[, c("rt", "mz", "Name")]

  #Candidate screening
  RT <- check_RT(feature_data = feature_data , reference =reference ) #check Retention Time similarity
  mass <- check_mass(feature_data = feature_data , reference = reference) #Calculating PPM error

  RT_mass <- check_RT_mass (RT , mass) # Merging ideal candidate in term of Retention Time and PPM error

  RI_sample <- RI_sample_data (feature_data_spectra=  feature_data_spectra) #establishing new sample data by calculating relative intensity
  RI_lib <- RI_lib_data ( reference_spectra =  reference_spectra , RT_mass)  ##establishing new library data by calculating relative intensity


  mass_mass <- check_mass_mass(RT_mass, RI_lib , RI_sample) #calculating and selecting proper PPM error and RI for mass mass


  return(list(candidate_list=mass_mass, RI_data=RI_sample, RI_reference=RI_lib))
}
