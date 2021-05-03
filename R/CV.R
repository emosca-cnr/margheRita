#' Select QC_samples and samples.
#' "QC" and "m" that indicates quality control samples and biological replicates in the title were used to select, respectively
#' @param df_post split into df_QC (dataframe of quality control) and df_samples (df of samples)
#' @import dplyr
#' @export



#apply CV QC vs CV of samples: eliminate when CV of QC >CV of samples for each metabolites

CV <- function(m_list, dirout) {
  dirout= paste(dirout, sep = "")
  dir.create(dirout)

  #subsetting metadata and df for QC and samples
  # we should discuss and agree on the followng code...
  m_list$QC_ann<- subset(m_list$sample_ann, m_list$sample_ann$type=="QC")
  m_list$QC<-m_list$df[colnames(m_list$df) %in% m_list$QC_ann$description]
  m_list$samples_ann<-subset(m_list$sample_ann, m_list$sample_ann$type=="sample")
  m_list$samples<-m_list$df[colnames(m_list$df) %in% m_list$samples_ann$description]
  CV_QC = c()
  for (i in 1:nrow(m_list$QC)) {
    mean_QC <- mean(as.numeric(m_list$QC[i, ])) #mean(col(QC))
    sd_QC <- sd(m_list$QC[i, ]) #sd(col(QC))
    CV_QC <- c(CV_QC, ((sd_QC / mean_QC) * 100))
    #CVQC= paste(dirout, sep = "")
    #utils::write.csv(CVQC, file = paste(dirout, "CVQC.csv"))
  }
  utils::write.csv(CV_QC, file = paste(dirout,"/CV_QC.csv", sep=""))
  CV_Samples <- c()
  for (j in 1:nrow(m_list$samples)) {
    mean_Samples <- mean(as.numeric(m_list$samples[j,])) #mean(col(Samples))
    sd_Samples <- sd(m_list$samples[j,]) #sd(col(Samples))
    CV_Samples <- c(CV_Samples, ((sd_Samples / mean_Samples) * 100))
    #CVSamples= paste(dirout, sep = "")
  }
  utils::write.csv(CV_Samples, file = paste(dirout,"/CV_Sample.csv", sep=""))
  CV_all <- cbind(CV_QC, CV_Samples)
  colnames(CV_all) = c("CV_QC", "CV_Samples")
  utils::write.csv(CV_all, paste(dirout, "/CV_all.csv", sep = ""))
  m_list$dfclean= m_list$dfclean[CV_QC < CV_Samples, ] #new m_list$df clean compared with the original

  return(m_list)
}


