#' Collapse
#' Subsetting metadata and data creating QC and samples data and metadata
#' Collapse technical replicates for samples.
#' @param For each sample mean was calculated to collapse technical replicates
#' @importFrom stats aggregate
#' @export


collapse<-function(m_list){
  m_list$sample_ann$pasted<-as.factor(paste(m_list$sample_ann$class,m_list$sample_ann$subclass,m_list$sample_ann$biological_rep,sep="_"))
   m_list$data<-as.data.frame(
    t(stats::aggregate(
      t(m_list$data), list(m_list$sample_ann$pasted), mean
    )))
  colnames(m_list$data) = m_list$data[1,]
 m_list$data <- m_list$data[-1, ]
   m_list$data = apply(m_list$data, MARGIN=1,as.numeric)
  m_list$data<-as.data.frame(t(m_list$data))

    return(m_list)
}
