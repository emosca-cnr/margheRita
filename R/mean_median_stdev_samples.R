#' Calculate mean, median and standard deviation for each biological group.
#'
#' @description Calculate mean, median and standard deviation for each biological group under study. This function works after collapse_tech_rep.
#'
#' @param mRList mRList object
#' @importFrom stats aggregate
#' @importFrom utils write.csv
#' @export
#' @return mRList object with mean, median and sd written in mRList$metab_ann element
#' @param dirout output directory

mean_median_stdev_samples<-function(mRList=NULL, dirout="./"){

  cat("According to dataset size, this might take a few minutes.\n")

  cat("Calculating means...\n")
  class_means <- t(apply(mRList$data, 1, function(x) tapply(x, mRList$sample_ann$class, mean)))
  colnames(class_means) <- paste0("mean_",  colnames(class_means))
  class_means <- data.frame(Feature_ID=rownames(class_means), class_means, stringsAsFactors = F)
  write.csv(class_means, file.path(dirout, "class_means.csv"), row.names = F)
  
  # mean_sample<-aggregate.data.frame(t(mRList$data), list(mRList$sample_ann$class), mean)
  # rownames(mean_sample)<-mean_sample$Group.1
  # mean_sample<-mean_sample[,-1]
  # mean_sample<-t(mean_sample)
  # mean_sample<-as.data.frame(mean_sample)
  # for (i in 1:length(names(mean_sample))) {
  #   names(mean_sample)[i] <- paste(names(mean_sample)[i],"mean", sep = "_")
  # }

  cat("Calculating medians...\n")
  class_medians <- t(apply(mRList$data, 1, function(x) tapply(x, mRList$sample_ann$class, median)))
  colnames(class_medians) <- paste0("median_",  colnames(class_medians))
  class_medians <- data.frame(Feature_ID=rownames(class_medians), class_medians, stringsAsFactors = F)
  write.csv(class_medians, file.path(dirout, "class_medians.csv"), row.names = F)
  
  # median_sample<-aggregate.data.frame(t(mRList$data), list(mRList$sample_ann$class), median)
  # rownames(median_sample)<-median_sample$Group.1
  # median_sample<-median_sample[,-1]
  # median_sample<-t(median_sample)
  # median_sample<-as.data.frame(median_sample)
  # for (i in 1:length(names(median_sample))) {
  #   names(median_sample)[i] <- paste(names(median_sample)[i],"median", sep = "_")
  # }

  cat("Calculating standard deviations...\n")
  class_sd <- t(apply(mRList$data, 1, function(x) tapply(x, mRList$sample_ann$class, sd)))
  colnames(class_sd) <- paste0("sd_",  colnames(class_sd))
  class_sd <- data.frame(Feature_ID=rownames(class_sd), class_sd, stringsAsFactors = F)
  write.csv(class_sd, file.path(dirout, "class_sd.csv"), row.names = F)
  
    # sd_sample<-aggregate.data.frame(t(mRList$data), list(mRList$sample_ann$class), sd)
  # rownames(sd_sample)<-sd_sample$Group
  # sd_sample<-sd_sample[,-1]
  # sd_sample<-t(sd_sample)
  # sd_sample<-as.data.frame(sd_sample)
  # for (i in 1:length(names(sd_sample))) {
  #   names(sd_sample)[i] <- paste(names(sd_sample)[i],"sd", sep = "_")
  # }

  # Mean_median_stdev <- cbind(mean_sample, median_sample, sd_sample)
  # mRList$metab_ann <- cbind(mRList$metab_ann, Mean_median_stdev)
  # 
  # names(mRList$metab_ann) <- gsub("1", "mean", names(mRList$metab_ann))
  # names(mRList$metab_ann) <- gsub("2", "median", names(mRList$metab_ann))

}
