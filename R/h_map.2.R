#' heatmap
#' @param mRList mRList
#' @export
#' @import ComplexHeatmap
#' @importFrom grDevices dev.off png
#' @importFrom graphics plot
#' @import viridis


h_map <- function(mRList, dirout="./", col_ann=NULL, col=NULL, scale_features=TRUE, features=NULL, top=500, ...){

  #dirout = paste(dirout, sep = "")
  dir.create(dirout, showWarnings = F)

  if(is.null(col_ann)){
    class_levels <- levels(as.factor(mRList$sample_ann$class))
    col_ann <- setNames(viridis::turbo(length(class_levels)), class_levels)
    
  }
  
  column_ha <- HeatmapAnnotation(class=mRList$sample_ann$class, col = list(class=col_ann))


  if(scale_features){
    data_scaled<-t(scale(t(mRList$data)))
  }
  

  if(is.null(features)){
    top_var <- apply(mRList$data, 1, var)
    top_var <- order(top_var, decreasing = T)[1:top]
  }
  
  if(is.null(col)){
    col <- viridis::viridis(7)
  }
  
  grDevices::png(paste0(dirout, "/Heatmap.png"), width= 180, height = 180, units = "mm", res = 300)
  hm <- ComplexHeatmap::Heatmap(as.matrix(data_scaled[top_var, ]), top_annotation = column_ha, col=col, ...)
  draw(hm)
  dev.off()

}
