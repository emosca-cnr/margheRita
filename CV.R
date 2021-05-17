#' Select QC_samples and samples.
#' "QC" and "m" that indicates quality control samples and biological replicates in the title were used to select, respectively
#' @param df_post split into df_QC (dataframe of quality control) and df_samples (df of samples)
#' @import dplyr
#' @export



#apply CV QC vs CV of samples: eliminate when CV of QC >CV of samples for each metabolites

CV <- function(m_list, dirout) {
  dirout= paste(dirout, sep = "")
  dir.create(dirout)

  CV_QC = c()
  mean_QC<-apply(m_list$QC,MARGIN=1, mean)
  sd_QC<-apply(m_list$QC,MARGIN = 1,sd)
 CV_QC <- c(CV_QC, ((sd_QC / mean_QC) * 100))
 utils::write.csv(CV_QC, file = paste(dirout,"/CV_QC.csv", sep=""))


  CV_Samples <- c()
    mean_Samples <- apply(m_list$data, MARGIN=1, mean)
    sd_Samples<-apply(m_list$data,MARGIN = 1,sd)
  CV_Samples <- c(CV_Samples, ((sd_Samples / mean_Samples) * 100))
   utils::write.csv(CV_Samples, file = paste(dirout,"/CV_Sample.csv", sep=""))

 CV_all <- cbind(CV_QC,CV_Samples)
  colnames(CV_all) = c("CV_QC", "CV_Samples")
  utils::write.csv(CV_all, paste(dirout, "/CV_all.csv", sep = ""))
  m_list$data<- m_list$data[CV_QC < CV_Samples, ] #m_list$data cleaned
  #clean metab (file containing information about metabolites) according to CV


  return(m_list)
}


