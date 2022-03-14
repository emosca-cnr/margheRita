#' generate the data set for annotation
#'
#'
#'
annotation_dataset <- function(column=c("HILIC", "LipC8", "pZIC", "RPLong", "RPShort"), mode=c("POS", "NEG"), accept_RI = 10){

## library

  lib_data <- margheRita_library(column = "HILIC", mode = "POS", accept_RI = 10)


## sample
sample_data <- as.data.frame(read_xlsx(system.file("extdata", "sample.xlsx", package = "margheRita")))
sample_data <- sample_data[, c(1, 2, 3, 7)]
colnames(sample_data) <- c("Feature_ID","rt", "mz", "MS_MS")
sample_data <- as.data.frame(sample_data[!is.na(sample_data$MS_MS), ])

sample_spectra <- sapply(sample_data$MS_MS, function(x) strsplit(x, " "))
sample_spectra <- lapply(sample_spectra, function(x) lapply(strsplit(x, ":"), as.numeric))
sample_spectra <- lapply(sample_spectra, function(x) do.call(rbind, x))
names(sample_spectra) <- sample_data$Feature_ID
sample_data$MS_MS <- NULL


return(list(lib_precursors =lib_data$lib_precursors, lib_peaks=lib_data$lib_peaks, lib_peaks_data=lib_data$lib_peaks_data, sample_data=sample_data, sample_spectra = sample_spectra))

}
