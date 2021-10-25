#' heatmap
#' @param m_list margheRita m_list
#' @export
#' @import ComplexHeatmap
#' @importFrom grDevices dev.off png
#' @importFrom graphics plot


h_map<-function(m_list, dirout="./"){

  #dirout = paste(dirout, sep = "")
  dir.create(dirout)

  column_ha <- HeatmapAnnotation(class= m_list$sample_ann$class, biorep= m_list$sample_ann$biological_rep,which="column")

  col_class <- as.factor(m_list$sample_ann$class)
  col_pal <- rainbow(length(levels(col_class)))
  col_biorep <- as.factor(m_list$sample_ann$biological_rep)
  col_bio <- rainbow(length(levels(col_biorep)))
  data_scaled<-t(scale(t(m_list$data)))
  #pdf(paste(dirout,"Heatmap.pdf", sep = ""))

  out_file <- paste0(dirout, "./Heatmap.png")
  grDevices::png(
    Heatmap,
    width= 8,
    height = 8,
    units = "in",
    res = 300
    )

  Heatmap(as.matrix(data_scaled[1:5, ]), top_annotation = column_ha,col=col_pal)

  grDevices::dev.off()

}
