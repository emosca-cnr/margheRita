#' mean, median, and stdev for technical replicates of each biological sample
#' @param metadata split in metadata for QC and metadata for samples
#' @param calculate mean, median and stdev for samples

#' @importFrom stats aggregate
#' @export


mean_media_stdev_samples<-function(m_list,dirout){
  dirout= paste(getwd(), "/mean_median/", sep = "")
  dir.create(dirout)

  #subsetting metadata and df for QC and samples the same in CV funtion

  #collapse only samples, we need to subset metadata only for samples
  #subsetting using column "type" in metadata QC" and "metadata_samples"
  #mean, median and stdev is calculated for QC samples and for samples taking into account only biological replicates
  #technical replicates of each biological replicates were collapsed

   #new column in metadata file
  m_list$samples_ann$pasted<-as.factor(paste(m_list$sample_ann$subclass,m_list$sample_ann$biological_rep,sep="_"))

  df_samples_mean <-
    (t(stats::aggregate(
      t(m_list$samples), list(m_list$samples_ann$pasted), mean
    )))
  colnames(df_samples_mean) = df_samples_mean[1, ]
  df_samples_mean <- df_samples_mean[-1, ]
  df_samples_mean = as.numeric(df_samples_mean)
  df_samples_mean = cbind(m_list$metabo_ann, df_samples_mean)

  df_samples_median <-
    (t(stats::aggregate(
      t(m_list$samples), list(m_list$samples_ann$pasted), median
    )))
  colnames(df_samples_median) = df_samples_median[1, ]
  df_samples_median <- df_samples_median[-1, ]
  df_samples_median = as.numeric(df_samples_median)
  df_samples_median = cbind(m_list$samples_ann, df_samples_median)

  df_samples_sd <-
    (t(stats::aggregate(
      t(m_list$samples), list(m_list$samples_ann$pasted), sd
    )))
  colnames(df_samples_sd) = df_samples_sd[1, ]
  df_samples_sd <- df_samples_sd[-1, ]
  df_samples_sd = as.numeric(df_samples_sd)
  df_samples_sd = cbind(m_list$samples_ann, df_samples_sd)


  Mean_median_stdev<-cbind( df_samples_mean,df_samples_median,df_samples_sd)
  media_median= paste(getwd(), "/mean_median/mean_median_stdev.csv", sep = "")
  utils::write.csv(Mean_median_stdev,media_median)
  return(list(m_list,df_samples_mean,df_samples_median,df_samples_sd))
  }
