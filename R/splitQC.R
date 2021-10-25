#' Split dataframe and metadata in two: one for samples and one for QC
#' @param m_list split metadata and dataframe for QC and for samples
#' @export
#' @return


splitQC<-function(m_list){

  idx_QC <- m_list$sample_ann$class=="QC"

  m_list$QC <- m_list$data[, idx_QC]
  m_list$QC_ann <- m_list$sample_ann[idx_QC, ]

  m_list$data <- m_list$data[, !idx_QC]
  m_list$sample_ann <- m_list$sample_ann[!idx_QC, ]


  return(m_list)

}

