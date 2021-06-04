
#'Collapse technical replicates for samples.
#'@param For each sample mean was calculated to collapse technical replicates
#' @importFrom stats aggregate
#' @export


collapse<-function(m_list){

  m_list$sample_ann$pasted<-as.factor(paste(m_list$sample_ann$class,m_list$sample_ann$subclass,m_list$sample_ann$biological_rep,sep="_"))
  m_list$data<-stats::aggregate.data.frame(t(m_list$data), list(m_list$sample_ann$pasted), mean)
  rownames(m_list$data)<-m_list$data$Group.1
  m_list$data<-m_list$data[,-1]
  m_list$data<-t(m_list$data)
  m_list$data<-as.data.frame(m_list$data)
  row.names(m_list$data)<-m_list$metab_ann$MS.Dial.ID

  target <- c(colnames(m_list$data))
  m_list$sample_ann<-m_list$sample_ann[match(target,m_list$sample_ann$pasted),]

  (m_list)
}
