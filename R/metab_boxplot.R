#' @importFrom graphics boxplot
#' @param m_list margheRita m_list
#' @export
#' @importFrom grDevices dev.off png



metab_boxplot<-function(m_list, dirout="./"){

  dirout = paste(dirout, sep = "")
  dir.create(dirout)

  col_factor <- as.factor(m_list$sample_ann$class)
  col_pal <- rainbow(length(levels(col_factor)))

  if (m_list$data$uni_corrected<0.05){
  data<-t(m_list$data)
  group<-as.factor(m_list$sample_ann$class)
  for (i in 1:ncol(data)) {
      boxplot(
     data[,i] ~ group , data=data,
      names =levels(group),
      #main = m_list$metabo_ann, #metabolites or Ms ID or annotation
      ylab="Relative Abundance", col= col_pal)
    png(file=paste(dirout,"metabolite",i,".png",sep=""), width = 480, height = 480, units = "px",
        bg = "white")
  }
  dev.off()
}
}
