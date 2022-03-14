#' a fake function to test pieces of code
#'
#'

margheRita_test_metab_annotation <- function(wdir="./"){

  ### ANNOTATION
  #data4annot <- generate_dataset_for_annotation()

  ### input data set
  input_data_file <- system.file("extdata", "dataset_drift.xlsx", package = "margheRita")
  input_metadata_file <- system.file("extdata", "dataset_drift_metadata.xlsx", package = "margheRita")
  mRList_raw <- read_input_file(input_data_file, metadata = input_metadata_file, data_start_col = 8, rt_col = 2, mz_col = 3)
  colnames(mRList_raw$metab_ann)[7] <- "MS_MS_spectrum"

  ### reference library
  #ref_lib <- as.data.frame(read_xlsx(system.file("extdata", "small_library.xlsx", package = "margheRita")))
  #ref_lib <- ref_lib[, c(2, 3, 4, 9, 32)]
  #colnames(ref_lib) <- c("rt", "mz", "Name", "rt_ref", "MS_MS_spectrum")
  #ref_lib <- ref_lib[ref_lib$Name != "Unknown", ]
  #ref_lib <- as.data.frame(ref_lib[!is.na(ref_lib$MS_MS_spectrum), ])

  ###LOAD the library
  lib_data <- select_library(column = "RPLong", mode = "POS", RI_min = 10)

  ###
  #ann_list <- metabolite_annotation(mRList = mRList_raw, reference = ref_lib)

  ###everyr step
  RT <-  check_RT(feature_data = mRList_raw$metab_ann , reference = lib_data$lib_precursor, rt_err_thr= 2)
  mass <- check_mass(feature_data = mRList_raw$metab_ann , reference = lib_data$lib_precursor, unaccept_flag=15, accept_flag=5, suffer_flag=10)
  RT_mass <- check_RT_mass (RT, mass, reference=lib_data$lib_precursor)

  #obrain spectra list from a string of mz intensity pairs separated by space, e.g. "55.01701:20 57.03321:163 57.8857:153 58.06513:9531"
  feature_spectra_list <- get_spectra_list_from_vector(spectra = setNames(mRList_raw$metab_ann$MS_MS_spectrum, mRList_raw$metab_ann$Feature_ID))
  #ref_spectra_list <- get_spectra_list_from_vector(spectra = setNames(ref_lib$MS_MS_spectrum, ref_lib$Name))

  ### calculation of Relative intensitities
  RI_sample <- RI_sample_data (feature_spectra=feature_spectra_list, acceptable_RI = 10 )
  #RI_lib <- RI_lib_data (reference = ref_spectra_list, RT_mass = RT_mass, acceptable_RI = 10)

  intense_peak <- check_intense_peak(RT_mass = RT_mass, RI_lib = lib_data$lib_peaks, reference = lib_data$lib_precursor, lib_peaks_cas=lib_data$lib_peaks_cas, mode = "POS", RI_sample = RI_sample, n_peaks=1, acceptable_PPM_err = 10)

  print(head(intense_peak))

}
