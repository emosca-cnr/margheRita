#' Subset
#' @param metadata split in metadata for QC and metadata for samples
#' @param calculate mean for each biological replicates
#' @importFrom stats aggregate
#' @export


subset<-function(m_list){
   m_list$QC_ann<- subset.data.frame(m_list$sample_ann, m_list$sample_ann$type=="QC")
  m_list$QC<-m_list$data[colnames(m_list$data) %in% m_list$QC_ann$description]
  m_list$sample_ann<-subset.data.frame(m_list$sample_ann, m_list$sample_ann$type=="sample")
   m_list$data<-m_list$data[colnames(m_list$data) %in% m_list$sample_ann$description]
   m_list$data<-as.data.frame(apply(m_list$data[2:nrow(m_list$data),],2,as.numeric))

  return(m_list)
}

