#' Metabolite Annotation
#'
#'
#' @export
#'
check_annotation = function(feature_data = NULL , reference = NULL , feature_data_spectra = NULL, reference_spectra= NULL){


  RT = check_RT(feature_data = feature_data , reference =reference ) #check Retention Time similarity
  mass = check_mass(feature_data = feature_data , reference = reference) #Calculating PPM error

  RT_mass = check_RT_mass (RT , mass) # Merging ideal candidate in term of Retention Time and PPM error

  RI_sample = RI_sample_data (feature_data_spectra=  feature_data_spectra) #establishing new sample data by calculating relative intensity
  RI_lib = RI_lib_data ( reference_spectra =  reference_spectra , RT_mass)  ##establishing new library data by calculating relative intensity


  mass_mass = check_mass_mass(RT_mass, RI_lib , RI_sample) #calculating and selecting proper PPM error and RI for mass mass


  return(mass_mass)
}
