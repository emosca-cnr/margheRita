#' margherRita - main function to run the full pipeline
#'
#'
#' @import MSnbase

generate_dataset_for_annotation <- function(wdir="./"){

  metabs_with_rt <- system.file("extdata", "serum_drift_GnpsMgf_0_2021371052.mgf", package = "margheRita")
  ref_lib_level2 <- system.file("extdata", "PSU-MSMLS.mgf", package = "margheRita")

  metabs_with_rt <- MSnbase::readMgfData(metabs_with_rt)
  ref_lib_level2 <- MSnbase::readMgfData(ref_lib_level2)

  #observed dataset with RT
  #print(head(fData(metabs_with_rt)))
  #cat(dim(fData(metabs_with_rt)))

  #metabs_with_rt@assayData$X1
  #print(as.data.frame(metabs_with_rt@assayData$X1))
  #mz(metabs_with_rt@assayData$X1)
  #intensity(metabs_with_rt@assayData$X1)

  #an example
  #plot(mz(metabs_with_rt@assayData$X1011), intensity(metabs_with_rt@assayData$X1011), type="h")
  #abline(v=precursorMz(metabs_with_rt@assayData$X1011), col="red")

  ### REFERENCE LIBRARY
  #print(head(fData(ref_lib_level2)))
  #cat(dim(fData(ref_lib_level2)))

  #ref_lib_level2@assayData$X1
  #print(as.data.frame(ref_lib_level2@assayData$X1))
  #mz(ref_lib_level2@assayData$X1)
  #intensity(ref_lib_level2@assayData$X1)

  ## REF LIBRARY "small library"
  ref_lib <- read_xlsx(system.file("extdata", "small_library.xlsx", package = "margheRita"))
  ref_lib <- ref_lib[, c(2, 3, 4, 9, 32)]
  colnames(ref_lib) <- c("rt", "mz", "Name", "rt_ref", "MS_MS")
  ref_lib <- ref_lib[ref_lib$Name != "Unknown", ]
  ref_lib <- as.data.frame(ref_lib[!is.na(ref_lib$MS_MS), ])
  ref_lib_spectra <- sapply(ref_lib$MS_MS, function(x) strsplit(x, " "))
  ref_lib_spectra <- lapply(ref_lib_spectra, function(x) lapply(strsplit(x, ":"), as.numeric))
  ref_lib_spectra <- lapply(ref_lib_spectra, function(x) do.call(rbind, x))
  names(ref_lib_spectra) <- ref_lib$Name
  ref_lib$MS_MS <- NULL

  #Sample data
  sample_data <- data.frame(Feature_ID=MSnbase::featureNames(metabs_with_rt), mz=as.numeric(fData(metabs_with_rt)[, "PEPMASS"]), rt=as.numeric(fData(metabs_with_rt)[, "RTINMINUTES"]), stringsAsFactors = F)

  sample_data_spectra <- lapply(MSnbase::spectra(metabs_with_rt), as.data.frame)
  sample_data_spectra <- sample_data_spectra[match(sample_data$ID, names(sample_data_spectra))]

  return(list(sample_data=sample_data, sample_data_spectra=sample_data_spectra, lib_data=ref_lib, lib_spectra=ref_lib_spectra))

}
