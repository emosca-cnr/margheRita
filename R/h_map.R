#' heatmap
## @color and names of samples depends on m_list$sample_ann$class
#' @param m_list margheRita m_list
#' @export
#' @import ComplexHeatmap
#' @importFrom grDevices dev.off pdf


h_map<-function(m_list, dirout="./"){

  #dirout = paste(dirout, sep = "")
  dir.create(dirout)
  column_ha <- HeatmapAnnotation(class= m_list$sample_ann$class, biorep= m_list$sample_ann$biological_rep)

  pdf(paste(dirout,"Heatmap.pdf", sep = ""))
  Heatmap(as.matrix(m_list$data[1:5, ]), top_annotation = column_ha)
  dev.off()

}
