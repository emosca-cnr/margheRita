#' Draw an heatmap
#' @param mRList mRList object
#' @param col_ann set of colors for sample annotation
#' @param column_ann sample_ann column for sample annotation
#' @param col set of color for data values 
#' @param scale_features whether to scale features or not
#' @param features names of features to plot (optiooal)
#' @param samples samples to consider (optional)
#' @param top only the top most variable features are plotted (if features is NULL)
#' @param ... further arguments for ComplexHeatmap::Heatmap
#' @export
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom grDevices dev.off png
#' @importFrom graphics plot
#' @importFrom stats var setNames
#' @importFrom pals brewer.rdylbu brewer.purples polychrome


h_map <- function(mRList=NULL, column_ann="class", col_ann=NULL, col=NULL, scale_features=TRUE, features=NULL, samples=NULL, top=20, ...){
  
  
  if(!is.null(samples)){
    data <- mRList$data[, colnames(mRList$data) %in% samples]
  }else{
    data <- mRList$data
  }

  
  if(scale_features){
    data<-t(scale(t(data)))
    hm_name <- "z"
  }else{
    hm_name <- "y"
  }
  
  if(is.null(features)){
    top_var <- apply(data, 1, var)
    top_var <- order(top_var, decreasing = T)[1:top]
    data <- data[top_var, ]
  }else{
    data <- data[rownames(data) %in% features, ]
    
  }
  
  
  column_ann <- as.factor(mRList$sample_ann[match(colnames(data), rownames(mRList$sample_ann)), column_ann])
  
  if(is.null(col_ann)){
    col_ann <- setNames(polychrome(length(levels(column_ann))), levels(column_ann))
  }
  
  column_ha <- HeatmapAnnotation(class=column_ann, col = list(class=col_ann))
  
  
  
  if(is.null(col)){
    if(scale_features){
      col <- rev(brewer.rdylbu(7))
    }else{
      col <- brewer.purples(7)
    }
  }
  
  Heatmap(as.matrix(data), top_annotation = column_ha, col=col, name = hm_name, ...)

}
