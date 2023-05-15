#' Calculate mean, median and standard deviation for each biological group.
#'
#' @description Calculate mean, median and standard deviation for each biological group under study. This function works after collapse_tech_rep.
#'
#' @param write_output (default FALSE), if TRUE table, as .csv format, is written and saved.
#' @param mRList mRList object
#' @importFrom stats aggregate
#' @importFrom utils write.csv
#' @export
#' @return mRList object with mean, median and sd written in mRList$metab_ann element
#' @param dirout output directory
#' @examples
#' ##library(dataset.margheRita)
#' ##dataset(norm_pos)
#' mRList<-mean_median_stdev_samples(mRList, write_output=TRUE)


mean_median_stdev_samples<-function(mRList, dirout, write_output=FALSE){

  cat("According to dataset size, this might take a few minutes.\n")


  mean_sample<-stats::aggregate.data.frame(t(mRList$data), list(mRList$sample_ann$class), mean)
  rownames(mean_sample)<-mean_sample$Group.1
  mean_sample<-mean_sample[,-1]
  mean_sample<-t(mean_sample)
  mean_sample<-as.data.frame(mean_sample)
  for (i in 1:length(names(mean_sample))) {
    names(mean_sample)[i] <- paste(names(mean_sample)[i],"mean", sep = "_")
  }

  median_sample<-stats::aggregate.data.frame(t(mRList$data), list(mRList$sample_ann$class), median)
  rownames(median_sample)<-median_sample$Group.1
  median_sample<-median_sample[,-1]
  median_sample<-t(median_sample)
  median_sample<-as.data.frame(median_sample)
  for (i in 1:length(names(median_sample))) {
    names(median_sample)[i] <- paste(names(median_sample)[i],"median", sep = "_")
  }


  sd_sample<-stats::aggregate.data.frame(t(mRList$data), list(mRList$sample_ann$class), sd)
  rownames(sd_sample)<-sd_sample$Group
  sd_sample<-sd_sample[,-1]
  sd_sample<-t(sd_sample)
  sd_sample<-as.data.frame(sd_sample)
  for (i in 1:length(names(sd_sample))) {
    names(sd_sample)[i] <- paste(names(sd_sample)[i],"sd", sep = "_")
  }



  Mean_median_stdev<-cbind(mean_sample,median_sample, sd_sample)
  mRList$metab_ann <- cbind(mRList$metab_ann, Mean_median_stdev)

  names(mRList$metab_ann) <- gsub("1", "mean", names(mRList$metab_ann))
  names(mRList$metab_ann) <- gsub("2", "median", names(mRList$metab_ann))


  if(write_output){
    dir.create(dirout)
    All= paste(dirout, "/mean_median_stdev.csv", sep ="")
    utils::write.csv(Mean_median_stdev, All)
    media= paste(dirout, "/mean.csv", sep ="")
    utils::write.csv(mean_sample, media)
    mediana= paste(dirout, "/median.csv", sep ="")
    utils::write.csv(median_sample, mediana)
    stdev= paste(dirout, "/stdev.csv", sep ="")
    utils::write.csv(sd_sample, stdev)
  }
  return(mRList)
}
