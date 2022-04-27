#' sample selection
#'
#' @export
#'
select_sample <- function(wdir="./"){


  ##batch1_norm_neg
  input_data_file <- as.data.frame(read_xlsx(system.file("extdata","dataset_drift.xlsx", package = "margheRita")))
  input_metadata_file <- as.data.frame(read_xlsx(system.file("extdata","dataset_drift_metadata.xlsx", package = "margheRita")))

  mRList_raw <- read_input_file(input_data_file, metadata = input_metadata_file, data_start_col = 7, rt_col = 2, mz_col = 3, MS_MS_column = 6)

  sample_data <- sample_data[, c(1, 2, 3, 5)]
  colnames(sample_data) <- c("Feature_ID","rt", "mz", "MS_MS")
  sample_data <- as.data.frame(sample_data[!is.na(sample_data$MS_MS), ])

  sample_spectra <- sapply(sample_data$MS_MS, function(x) strsplit(x, " "))
  sample_spectra <- lapply(sample_spectra, function(x) lapply(strsplit(x, ":"), as.numeric))
  sample_spectra <- lapply(sample_spectra, function(x) do.call(rbind, x))
  names(sample_spectra) <- sample_data$Feature_ID
  sample_data$MS_MS <- NULL

  return(list(feature_data= sample_data, feature_spectra = sample_spectra))
}
