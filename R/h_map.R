#' Draw an heatmap
#' @param mRList mRList object
#' @param col_ann set of colors for sample annotation
#' @param data.use whether to use mRList$data (default) or mRList$data_ann, which is the annotated version obtained via metabolite_identification
#' @param column_ann sample_ann column for sample annotation
#' @param col set of color for data values
#' @param scale_features whether to scale features or not
#' @param features names of features to plot (optional)
#' @param feature_label "Feature_ID" or "Name" from the mRList$data_ann
#' @param samples samples to consider (optional)
#' @param top only the top most variable features are plotted (if features is NULL)
#' @param use_top_annotated_metab if TRUE, it'll use the best name for the metabolite from the annotation, considering the level and also the top significant specified in the following arguments.
#' @param test_method the name of the statistical test table contained in the mRlist object
#' @param test_value the column of the statistical test table to consider (e.g.: p or q)
#' @param ... further arguments for ComplexHeatmap::Heatmap
#' @export
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom grDevices dev.off png
#' @importFrom graphics plot
#' @importFrom stats var setNames
#' @importFrom pals brewer.rdylbu brewer.purples polychrome



h_map <- function(mRList=NULL, column_ann="class", data.use = c("data", "data_ann"), col_ann=NULL, col=NULL, scale_features=TRUE, features=NULL, feature_label = c("Feature_ID", "Name"), samples=NULL, top=200, use_top_annotated_metab = FALSE, test_method="anova", test_value = "q", ...){
  
  data.use <- match.arg(data.use, c("data", "data_ann"))
  test_value <- match.arg(test_value, c("q", "p"))
  feature_label <- match.arg(feature_label, c("Feature_ID", "Name"))
  
  data <- mRList[[data.use]]
  stopifnot(length(data)>0)
  
  # set rownames and remove character columns
  if(data.use == "data_ann"){
    
    if (feature_label == "Feature_ID") {
      rownames(data) <- data$Feature_ID
    }
    if (feature_label == "Name") {
      rownames(data) <- data$Name
    } 
    
    data <- data[, ! colnames(data) %in% c("Feature_ID", "Name")]
  }
  
  #Selecting the best annotation metabolite as a name, also considering the statistics:
  if (use_top_annotated_metab) {
    
    if (data.use != "data_ann") {
      stop('use_top_annotated_metab=TRUE requires data.use="data_ann"')
    }
    
    top_sig_feat <- rownames(mRList[[test_method]])[order(mRList[[test_method]][, test_value])] #get results and order by test_value
    if (feature_label == "Name") { #Name instead of Feature_ID
      top_sig_feat <- mRList$data_ann$Name[match(top_sig_feat, mRList$data_ann$Feature_ID)]
      top_sig_feat <- top_sig_feat[!is.na(top_sig_feat)]
    }
    
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
    }else{
      col <- brewer.purples(7)
    }
  }
  
  print(Heatmap(as.matrix(data), top_annotation = column_ha, col=col, name = hm_name, ...))
  
}

