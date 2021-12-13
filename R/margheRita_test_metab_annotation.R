#' margherRita - main function to run the full pipeline
#'
#'
#' @import notame MSnbase

margheRita_test_metab_annotation <- function(wdir="./"){

  ### ANNOTATION
  data4annot <- generate_dataset_for_annotation()

  ### input data set
  input_data_file <- system.file("extdata", "dataset_drift.xlsx", package = "margheRita")
  input_metadata_file <- system.file("extdata", "dataset_drift_metadata.xlsx", package = "margheRita")
  mRList_raw <- read_input_file(input_data_file, metadata = input_metadata_file, data_start_col = 8, rt_col = 2, mz_col = 3)
  colnames(mRList_raw$metab_ann)[7] <- "MS_MS_spectrum"

  ### reference library
  ref_lib <- as.data.frame(read_xlsx(system.file("extdata", "small_library.xlsx", package = "margheRita")))
  ref_lib <- ref_lib[, c(2, 3, 4, 9, 32)]
  colnames(ref_lib) <- c("rt", "mz", "Name", "rt_ref", "MS_MS_spectrum")
  ref_lib <- ref_lib[ref_lib$Name != "Unknown", ]
  ref_lib <- as.data.frame(ref_lib[!is.na(ref_lib$MS_MS_spectrum), ])


  ann_list <- metabolite_annotation(mRList = mRList_raw, reference = ref_lib)

}
