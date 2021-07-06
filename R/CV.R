#' Select QC_samples and samples.
#' "QC" and "m" that indicates quality control samples and biological replicates in the title were used to select, respectively
#' @param df_post split into df_QC (dataframe of quality control) and df_samples (df of samples)
#' @export



#apply CV QC vs CV of samples: eliminate when CV of QC >CV of samples for each metabolites

CV <- function(m_list, dirout) {

  dir.create(dirout)

  mean_QC <- apply(m_list$QC, MARGIN=1, mean)
  sd_QC <- apply(m_list$QC, MARGIN = 1, sd)
  CV_QC <- (sd_QC / mean_QC) * 100
  utils::write.csv(CV_QC, file = paste0(dirout,"/CV_QC.csv"))


  mean_Samples <- apply(m_list$data, MARGIN=1, mean)
  sd_Samples <- apply(m_list$data,MARGIN = 1,sd)
  CV_Samples <- (sd_Samples / mean_Samples) * 100
  utils::write.csv(CV_Samples, file = paste0(dirout,"/CV_Sample.csv"))

  CV_all <- cbind(CV_QC, CV_Samples)
  colnames(CV_all) = c("CV_QC", "CV_Samples")
  utils::write.csv(CV_all, paste0(dirout, "/CV_all.csv"))

  idx_keep <- CV_QC < CV_Samples
  cat("Metabolites with appropriate CV\n")
  print(table(idx_keep))

  m_list$data <- m_list$data[idx_keep, ] #m_list$data cleaned
  m_list$metab_ann <- m_list$metab_ann[idx_keep, ] #m_list$data cleaned

  m_list$QC <- m_list$QC[idx_keep, ] #m_list$data cleaned

  return(m_list)
}


