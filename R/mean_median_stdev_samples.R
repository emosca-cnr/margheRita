#' Calculate mean, median and standard deviation
#' Calculate mean, median and standard deviation for each level of class. It provides the results in text files in dirout directory.
#' @param mRList mRList object
#' @param dirout output directory
#' @importFrom stats aggregate median sd
#' @importFrom utils write.csv
#' @returns NULL
#' @export

mean_median_stdev_samples <- function(mRList=NULL, dirout=NULL){

  if (!is.null(dirout)) {
    dir.create(dirout, showWarnings = F)
  }else{
    dirout <- getwd()
  }
  
  cat("According to dataset size, this might take a few minutes.\n")

  cat("Calculating means...\n")
  class_means <- t(apply(mRList$data, 1, function(x) tapply(x, mRList$sample_ann$class, mean)))
  colnames(class_means) <- paste0("mean_",  colnames(class_means))
  class_means <- data.frame(Feature_ID=rownames(class_means), class_means, stringsAsFactors = F)
  write.csv(class_means, file.path(dirout, "class_means.csv"), row.names = F)

  cat("Calculating medians...\n")
  class_medians <- t(apply(mRList$data, 1, function(x) tapply(x, mRList$sample_ann$class, median)))
  colnames(class_medians) <- paste0("median_",  colnames(class_medians))
  class_medians <- data.frame(Feature_ID=rownames(class_medians), class_medians, stringsAsFactors = F)
  write.csv(class_medians, file.path(dirout, "class_medians.csv"), row.names = F)

  cat("Calculating standard deviations...\n")
  class_sd <- t(apply(mRList$data, 1, function(x) tapply(x, mRList$sample_ann$class, sd)))
  colnames(class_sd) <- paste0("sd_",  colnames(class_sd))
  class_sd <- data.frame(Feature_ID=rownames(class_sd), class_sd, stringsAsFactors = F)
  write.csv(class_sd, file.path(dirout, "class_sd.csv"), row.names = F)

}
