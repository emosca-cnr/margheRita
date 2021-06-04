#' Split the QC samples
#' @param metadata split in metadata for QC and metadata for samples
#' @param calculate mean for each biological replicates
#' @importFrom stats aggregate
#' @export


splitQC<-function(m_list){

  #m_list$QC_ann <- subset.data.frame(m_list$sample_ann, m_list$sample_ann$class=="QC") #ettore type->class
  #m_list$QC <- subset.data.frame(m_list$data, m_list$sample_ann$class=="QC") #ettore
  #m_list$QC < m_list$data[colnames(m_list$data) %in% m_list$QC_ann$description]
  #m_list$sample_ann <- subset.data.frame(m_list$sample_ann, m_list$sample_ann$class=="sample")
  #m_list$data <- m_list$data[colnames(m_list$data) %in% m_list$sample_ann$id]
  #m_list$data <- as.data.frame(apply(m_list$data[2:nrow(m_list$data),], 2, as.numeric))

  idx_QC <- m_list$sample_ann$class=="QC"

  m_list$QC <- m_list$data[, idx_QC]
  m_list$QC_ann <- m_list$sample_ann[idx_QC, ]

  m_list$data <- m_list$data[, !idx_QC]
  m_list$sample_ann <- m_list$sample_ann[!idx_QC, ]


  return(m_list)

}

