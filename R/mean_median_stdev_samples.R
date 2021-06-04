#' mean, median, and stdev for technical replicates of each biological sample
#' @param metadata split in metadata for QC and metadata for samples
#' @param calculate mean, median and stdev for samples

#' @importFrom stats aggregate, sd, mean, median
#' @export


mean_media_stdev_samples<-function(m_list,dirout){
  dirout = paste(dirout, sep = "")
  dir.create(dirout)

mean_sample<-stats::aggregate.data.frame(t(m_list$data), list(m_list$sample_ann$pasted), mean)
rownames(mean_sample)<-mean_sample$Group.1
mean_sample<-mean_sample[,-1]
mean_sample<-t(mean_sample)
mean_sample<-as.data.frame(mean_sample)

median_sample<-stats::aggregate.data.frame(t(m_list$data), list(m_list$sample_ann$pasted), median)
rownames(median_sample)<-median_sample$Group.1
median_sample<-median_sample[,-1]
median_sample<-t(median_sample)
median_sample<-as.data.frame(median_sample)


sd_sample<-stats::aggregate.data.frame(t(m_list$data), list(m_list$sample_ann$pasted), sd)
rownames(sd_sample)<-sd_sample$Group.1
sd_sample<-sd_sample[,-1]
sd_sample<-t(sd_sample)
sd_sample<-as.data.frame(sd_sample)



  Mean_median_stdev<-cbind(mean_sample,median_sample, sd_sample )
  media_median= paste(getwd(), "/mean_median/mean_median_stdev.csv", sep = "")
  utils::write.csv(Mean_median_stdev,media_median)
  return(m_list)
  }
