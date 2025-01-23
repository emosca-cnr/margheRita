#' Draw an heatmap
#' @param mRList mRList object
#' @param col_ann set of colors for sample annotation
#' @param data.use whether to use mRList$data (default) or mRList$data_ann, which is the annotated version obtained via metabolite_identification
#' @param column_ann sample_ann column for sample annotation
#' @param col set of color for data values
#' @param scale_features whether to scale features or not
#' @param features names of features to plot (optional)
#' @param feature_id "Feature_ID" or "Name" from the mRList$data_ann
#' @param samples samples to consider (optional)
#' @param top only the top most variable features are plotted (if features is NULL)
#' @param use_top_annoteted_metab if TRUE, it'll use the best name for the metabolite from the annotation, considering the level and also the top significant specified in the following arguments.
#' @param test_method the name of the statistical test table contained in the mRlist object
#' @param test_value the column of the statistical test table to consider (e.g.: p or q)
#' @param ... further arguments for ComplexHeatmap::Heatmap
#' @export
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom grDevices dev.off png
#' @importFrom graphics plot
#' @importFrom stats var setNames
#' @importFrom pals brewer.rdylbu brewer.purples polychrome



h_map <- function(mRList=NULL, column_ann="class", data.use = c("data", "data_ann"), col_ann=NULL, col=NULL, scale_features=TRUE, features=NULL, feature_id = "Feature_ID", samples=NULL, top=200, use_top_annoteted_metab = FALSE, test_method="anova", test_value = "q", ...){
  
  data.use <- match.arg(data.use, c("data", "data_ann"))
  
  data <- mRList[[data.use]]
  stopifnot(length(data)>0)
  if(data.use == "data_ann"){
    
    if (feature_id == "Feature_ID") {
      rownames(data) <- data$Feature_ID
    } else if (feature_id == "Name") {
      rownames(data) <- data$Name
    } else {
      stop('feature_id must be "Feature_ID" or "Name"')
    }
    
    data <- data[, ! colnames(data) %in% c("Feature_ID", "Name")]
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
  
  if(is.null(features)){
    if(length(top_var)<top) {top <- length(top_var)}
    top_var <- apply(data, 1, var)
    top_var <- order(top_var, decreasing = T)[1:top]
    data <- data[top_var, ]
  }else{
    if(length(features)>top) {features <- features[1:top]}
    data <- data[rownames(data) %in% features, ]
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
  
  #GF, selecting the best annotation metabolite as a name, also considering the statistics:
  
  if (use_top_annoteted_metab) {
    if (data.use == "data_ann") {
      
      new_row_names_data <- row.names(data)
      
      statistical_test_table <- mRList[[test_method]][order(as.vector(mRList[[test_method]][,test_value])),]
      statistical_test_table[,"Feature_ID"] <- row.names(statistical_test_table)
      statistical_test_table <- statistical_test_table[which(!duplicated(statistical_test_table$Feature_ID)),]
      
      annotation_table <- mL_norm[["metabolite_identification"]][["associations"]][order(as.numeric(as.vector(mL_norm[["metabolite_identification"]][["associations"]][,"Level"]))),]
      annotation_table <- annotation_table[which(!duplicated(annotation_table$Feature_ID)),]
      
      combined_table <- merge(x = statistical_test_table, y = annotation_table, by = "Feature_ID")
      
      
      
      for (i in 1:length(new_row_names_data)) {
        this_feat_id <- as.vector(mRList[[data.use]][,"Feature_ID"])[which(as.vector(mRList[[data.use]][,"Name"]) ==  row.names(data)[i])]
        
        if (length(this_feat_id) ==1) {
          if (this_feat_id %in% combined_table$Feature_ID) {
            new_row_names_data[i] <- combined_table$Name[which(combined_table$Feature_ID==this_feat_id)]
          }
        } else if (length(this_feat_id) >1) {
          this_feat_id <- this_feat_id[which(this_feat_id %in% combined_table$Feature_ID)]
          if (length(this_feat_id)>0) {
            new_row_names_data[i] <- combined_table$Name[which(combined_table$Feature_ID==this_feat_id)][1]
          }
        }
      }
      
      row.names(data) <- new_row_names_data
    }
  }
  
  
  row.names(data) <- sub(";.*", "", row.names(data))
  
  #png(filename = file.path(dirout, "Heatmap.png"), width = 200, height = 200, units = "mm", res=300)
  
  print(Heatmap(as.matrix(data), top_annotation = column_ha, col=col, name = hm_name, ...))
  
  #dev.off()
  
}

