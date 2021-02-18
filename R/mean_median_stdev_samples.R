#mean, median, and stdev for technical replicates of each biological sample
#@param metadata split in metadata for QC and metadata for samples
#@param calculate mean, median and stdev for samples

library(dplyr)
df_QC <- select(dataframe_post, contains("QC"))
df_samples <- select(dataframe_post, contains("m")) #

#collapse only samples, we need to subset metadata only for samples
#subsetting using column "type" in metadata QC" and "metadata_samples"
#mean, median and stdev is calculated for QC samples and for samples taking into account only biological replicates
#technical replicates of each biological replicates were collapsed

metadata_QC <- metadata[metadata$type == "QC", ]
metadata_samples <- metadata[metadata$type == "sample", ]

#new column in metadata file
metadata$pasted<-as.factor(paste(metadata$subclass,metadata$biological_rep,sep="_"))

df_samples_mean <-
  (t(aggregate(
    t(df_samples[, 4:ncol(df_samples)]), list(metadata_samples$pasted), mean
  )))
colnames(df_samples_mean) = df_samples_mean[1, ]
df_samples_mean <- df_samples_mean[-1, ]
df_samples_mean = as.numeric(df_samples_mean)
df_samples_mean = cbind(dataframe[, 1:3], df_samples_mean)

df_samples_median <-
  (t(aggregate(
    t(df_samples[, 4:ncol(df_samples)]), list(metadata_samples$pasted), median
  )))
colnames(df_samples_median) = df_samples_median[1, ]
df_samples_median <- df_samples_median[-1, ]
df_samples_median = as.numeric(df_samples_median)
df_samples_median = cbind(dataframe[, 1:3], df_samples_median)

df_samples_sd <-
  (t(aggregate(
    t(df_samples[, 4:ncol(df_samples)]), list(metadata_samples$pasted), sd
  )))
colnames(df_samples_sd) = df_samples_sd[1, ]
df_samples_sd <- df_samples_sd[-1, ]
df_samples_sd = as.numeric(df_samples_sd)
df_samples_sd = cbind(dataframe[, 1:3], df_samples_sd)
