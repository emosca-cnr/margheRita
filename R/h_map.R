#' Draw an heatmap
#' @param mRList mRList object
#' @param col_ann set of colors for sample annotation
#' @param data.use whether to use mRList$data (default) or mRList$data_ann, which is the annotated version obtained via metabolite_identification
#' @param column_ann sample_ann column for sample annotation
#' @param col set of color for data values
#' @param scale_features whether to scale features or not
#' @param features names of features to plot (optional)
#' @param samples samples to consider (optional)
#' @param top only the top most variable features are plotted (if features is NULL)
#' @param use_top_annotated_metab if TRUE, it'll use the best name for the metabolite from the annotation, considering the level and also the top significant specified in the following arguments.
#' @param test_method the name of the statistical test table contained in the mRlist object
#' @param test_value the column of the statistical test table to consider (e.g.: p or q)
#' @param truncate_row_labels FALSE or a integer that specifies the number of characters to keep
#' @param rm_outliers remove outliers when setting colours
#' @param ... further arguments for ComplexHeatmap::Heatmap
#' @export
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom grDevices dev.off png boxplot.stats
#' @importFrom graphics plot
#' @importFrom stats var setNames
#' @importFrom pals brewer.rdylbu brewer.purples polychrome
#' @importFrom stringr str_trunc
#' @importFrom circlize colorRamp2


h_map <- function(mRList=NULL, column_ann="class", data.use = c("data", "data_ann"), col_ann=NULL, col=NULL, scale_features=TRUE, features=NULL, samples=NULL, top=200, use_top_annotated_metab = FALSE, test_method="anova", test_value = "q", truncate_row_labels=FALSE, rm_outliers=FALSE, ...){
  
  data.use <- match.arg(data.use, c("data", "data_ann"))
  test_value <- match.arg(test_value, c("q", "p"))
  
  data <- mRList[[data.use]]
  stopifnot(length(data)>0)
  
  # set rownames and remove character columns
  if(data.use == "data_ann"){
    
    rownames(data) <- data$Name
    data <- data[, ! colnames(data) %in% c("Feature_ID", "Name")]
    
  }
  
  #Selecting the best annotation metabolite as a name, also considering the statistics:
  if (use_top_annotated_metab) {
    
    if (data.use != "data_ann") {
      stop('use_top_annotated_metab=TRUE requires data.use="data_ann"')
    }
    
    top_sig_feat <- rownames(mRList[[test_method]])[order(mRList[[test_method]][, test_value])] #get results and order by test_value
 
    data <- data[match(top_sig_feat, rownames(data)), ]
    
  }else{
    
    if(is.null(features)){ #by var
      top_var <- apply(data, 1, var)
      top_var <- order(top_var, decreasing = T)[1:max(top, length(top_var))]
      data <- data[top_var, ]
    }else{ #by given name
      data <- data[match(features[1:min(top, length(features))], rownames(data)), ]
    }
    
  }
  
  if(!is.null(samples)){
    data <- data[, colnames(data) %in% samples]
  }
  stopifnot(length(data)>0)
  
  if(scale_features){
    data<-t(scale(t(data)))
    hm_name <- "z"
  }else{
    hm_name <- "y"
  }
  
  
  column_ann <- as.factor(mRList$sample_ann[match(colnames(data), rownames(mRList$sample_ann)), column_ann])
  
  if(is.null(col_ann)){
    col_ann <- setNames(polychrome(length(levels(column_ann))), levels(column_ann))
  } else {
    if (length(col_ann) != length(levels(column_ann))) {stop("the length of col_ann must be identical to the length of groups")}
    names(col_ann) <- levels(column_ann)
  }
  
  column_ha <- HeatmapAnnotation(class=column_ann, col = list(class=col_ann))
  
  
  if(is.null(col)){
    if(scale_features){
      col <- rev(brewer.rdylbu(7))
      if(rm_outliers){
        min_max <- max(abs(boxplot.stats(as.numeric(data))$stats[c(1, 5)]))
        min_max <- c(-min_max, 0, min_max)
      }else{
        min_max <- c(min(data), 0, max(data))
      }
      col <- colorRamp2(min_max, col[c(1, 4, 7)])
      
    }else{
      col <- brewer.purples(7)
    }
  }
  
  if(is.numeric(truncate_row_labels)){
    row_labels <- str_trunc(rownames(data), truncate_row_labels)
  }else{
    row_labels <- rownames(data)
  }
  print(Heatmap(as.matrix(data), top_annotation = column_ha, col=col, name = hm_name, row_labels=row_labels, ...))
  
}

