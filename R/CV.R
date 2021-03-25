#' Select QC_samples and samples.
#' "QC" and "m" that indicates quality control samples and biological replicates in the title were used to select, respectively
#' @param df_post split into df_QC (dataframe of quality control) and df_samples (df of samples)
#' @import dplyr
#' @export


#create directory CV QC and CV samples
#CVdir = paste(getwd(), "/CV/", sep = "")
#dir.create(CVdir)

#apply CV QC vs CV of samples: eliminate when CV of QC >CV of samples for each metabolites

CV <- function(df_post) {
  #CVdir = paste(getwd(), "/CV/", sep = "")
  #dirout(CVdir)
  #create two matrices:one for QC and one for samples, in function?
  df_QC <- select(df_post, contains("QC"))
  df_samples <- select(df_post, contains("m"))
  CV_QC = c()
  for (i in 1:nrow(df_QC[,4:ncol(df_QC)])) {
    mean_QC <- mean(as.numeric(df_QC[i, ])) #mean(col(QC))
    sd_QC <- sd(df_QC[i, ]) #sd(col(QC))
    CV_QC <- c(CV_QC, ((sd_QC / mean_QC) * 100))
  }
  CV_Samples <- c()
  for (j in 1:nrow(df_samples)) {
    mean_Samples <- mean(as.numeric(df_samples[j,4:ncol(df_samples)])) #mean(col(Samples))
    sd_Samples <- sd(df_samples[j, ]) #sd(col(Samples))
    CV_Samples <- c(CV_Samples, ((sd_Samples / mean_Samples) * 100))
  }

  CV_all <- cbind(CV_QC, CV_Samples)
  colnames(CV_all) = c("CV_QC", "CV_Samples")
  write.csv(CV_all, paste(CVdir, "CV_all.csv", sep = ""))
  dataframe_post = dataframe_post[CV_QC < CV_Samples, ]
  return(dataframe_post)
}


