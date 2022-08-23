#' Draw an heatmap
#' @param mRList mRList object
#' @param dirout output directory
#' @param col_ann set of colors for sample annotation
#' @param col set of color for data values 
#' @param scale_features whether to scale features or not
#' @param features names of features to plot (optiooal)
#' @param samples samples to consider (optional)
#' @param top only the top most variable features are plotted (if features is NULL)
#' @param ... further arguments for ComplexHeatmap::Heatmap
#' @export
#' @import ComplexHeatmap
#' @importFrom grDevices dev.off png
#' @importFrom graphics plot
#' @import viridis


h_map <- function(mRList, dirout="./", col_ann=NULL, col=NULL, scale_features=TRUE, features=NULL, samples=NULL, top=500, ...){

  #dirout = paste(dirout, sep = "")
  dir.create(dirout, showWarnings = F)

  if(scale_features){
    data<-t(scale(t(mRList$data)))
  }
  

  if(is.null(features)){
    top_var <- apply(mRList$data, 1, var)
    top_var <- order(top_var, decreasing = T)[1:top]
    data <- data[top_var, ]
  }else{
    data <- data[rownames(data) %in% features, ]
    
  }
  
  if(!is.null(samples)){
    data <- data[, colnames(data) %in% samples]
  }

  
  if(is.null(col_ann)){
    class <- as.factor(mRList$sample_ann$class[rownames(mRList$sample_ann) %in% samples])
    col_ann <- setNames(viridis::turbo(length(levels(class))), levels(class))
    
  }
  
  column_ha <- HeatmapAnnotation(class=class, col = list(class=col_ann))
  
  
    
  if(is.null(col)){
    col <- viridis::viridis(7)
  }
  
  grDevices::png(paste0(dirout, "/Heatmap.png"), width= 180, height = 180, units = "mm", res = 300)
  hm <- ComplexHeatmap::Heatmap(as.matrix(data), top_annotation = column_ha, col=col, ...)
  draw(hm)
  dev.off()

}
