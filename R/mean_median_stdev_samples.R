#' mean, median, and stdev for technical replicates of each biological sample
#' @param metadata split in metadata for QC and metadata for samples
#' @param calculate mean, median and stdev for samples
#' @importFrom stats aggregate
#' @importFrom utils write.csv
#' @export


mean_media_stdev_samples<-function(m_list,dirout){

  cat("According to dataset size, this might take a few minutes.\n")
  dirout = paste(dirout, sep = "")
  dir.create(dirout)

  mean_sample<-stats::aggregate.data.frame(t(m_list$data), list(m_list$sample_ann$class), mean)
  rownames(mean_sample)<-mean_sample$Group.1
  mean_sample<-mean_sample[,-1]
  mean_sample<-t(mean_sample)
  mean_sample<-as.data.frame(mean_sample)
  for (i in 1:length(names(mean_sample))) {
     names(mean_sample)[i] <- paste(names(mean_sample)[i],"mean", sep = "_") 
  }

  median_sample<-stats::aggregate.data.frame(t(m_list$data), list(m_list$sample_ann$class), median)
  rownames(median_sample)<-median_sample$Group.1
  median_sample<-median_sample[,-1]
  median_sample<-t(median_sample)
  median_sample<-as.data.frame(median_sample)
  for (i in 1:length(names(median_sample))) {
     names(median_sample)[i] <- paste(names(median_sample)[i],"median", sep = "_") 
  }


  sd_sample<-stats::aggregate.data.frame(t(m_list$data), list(m_list$sample_ann$class), sd)
  rownames(sd_sample)<-sd_sample$Group
  sd_sample<-sd_sample[,-1]
  sd_sample<-t(sd_sample)
  sd_sample<-as.data.frame(sd_sample)
  for (i in 1:length(names(sd_sample))) {
     names(sd_sample)[i] <- paste(names(sd_sample)[i],"sd", sep = "_") 
  }



  Mean_median_stdev<-cbind(mean_sample,median_sample, sd_sample)
  m_list$metab_ann <- cbind(m_list$metab_ann, Mean_median_stdev)

  names(m_list$metab_ann) <- gsub("1", "mean", names(m_list$metab_ann))
  names(m_list$metab_ann) <- gsub("2", "median", names(m_list$metab_ann))


  All= paste(dirout, "/mean_median_stdev.csv", sep ="")
  utils::write.csv(Mean_median_stdev, All)
  media= paste(dirout, "/mean.csv", sep ="")
  utils::write.csv(mean_sample, media)
  mediana= paste(dirout, "/median.csv", sep ="")
  utils::write.csv(median_sample, mediana)
  stdev= paste(dirout, "/stdev.csv", sep ="")
  utils::write.csv(sd_sample, stdev)

  return(m_list)
}
