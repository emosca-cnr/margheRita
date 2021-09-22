#' margherRita - main function to run the full pipeline
#'
#'
#' @import notame MSnbase

generate_dataset_for_annotation <- function(wdir="./"){

  metabs_with_rt <- system.file("extdata", "serum_drift_GnpsMgf_0_2021371052.mgf", package = "margheRita")
  ref_lib_level2 <- system.file("extdata", "PSU-MSMLS.mgf", package = "margheRita")

  metabs_with_rt <- MSnbase::readMgfData(metabs_with_rt)
  ref_lib_level2 <- MSnbase::readMgfData(ref_lib_level2)

  #observed dataset with RT
  print(head(fData(metabs_with_rt)))
  cat(dim(fData(metabs_with_rt)))

  metabs_with_rt@assayData$X1
  print(as.data.frame(metabs_with_rt@assayData$X1))
  mz(metabs_with_rt@assayData$X1)
  intensity(metabs_with_rt@assayData$X1)

  #an example
  plot(mz(metabs_with_rt@assayData$X1011), intensity(metabs_with_rt@assayData$X1011), type="h")
  abline(v=precursorMz(metabs_with_rt@assayData$X1011), col="red")

  ### REFERENCE LIBRARY
  print(head(fData(ref_lib_level2)))
  cat(dim(fData(ref_lib_level2)))

  ref_lib_level2@assayData$X1
  print(as.data.frame(ref_lib_level2@assayData$X1))
  mz(ref_lib_level2@assayData$X1)
  intensity(ref_lib_level2@assayData$X1)


}
