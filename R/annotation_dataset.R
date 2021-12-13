#' generate the data set for annotation 
#'
#'
#' 
annotation_dataset <- function(wdir="./"){


## library
library_data <- as.data.frame(read_xlsx(system.file("extdata", "small_library.xlsx", package = "margheRita")))
library_data <- library_data[, c(2, 3, 4, 32)]
colnames(library_data) <- c("rt", "mz", "Name", "MS_MS")
library_data <- library_data[library_data$Name != "Unknown", ]
library_data <- as.data.frame(library_data[!is.na(library_data$MS_MS), ])

library_spectra <- sapply(library_data$MS_MS, function(x) strsplit(x, " "))
library_spectra <- lapply(library_spectra, function(x) lapply(strsplit(x, ":"), as.numeric))
library_spectra <- lapply(library_spectra, function(x) do.call(rbind, x))
names(library_spectra) <- library_data$Name
library_data$MS_MS <- NULL

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


return(list(lib_data=library_data, lib_spectra= library_spectra, sample_data=sample_data, sample_spectra = sample_spectra))

}