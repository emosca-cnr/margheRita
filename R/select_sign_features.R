#' select significant features
#' 
#' The output of this function can be passed in the features argument of the h_map function
#' 
#' @param mRList mRList object
#' @param test_method the name of the statistical test table contained in the mRlist object
#' @param test_value the column of the statistical test table to consider (e.g.: p or q)
#' @param cutoff_value the value of p-value cut-off to consider a feature as significant
#' @param feature_id "Feature_ID" or "Name" from the mRList$data_ann
#' @param values if TRUE a named vector with test_value values will be given in output; the beast value will be considered in case of 1-to-n mappings
#' @param split_by_semicolon if TRUE, the names will be split by semicolons.
#' @importFrom stats setNames
#' @return a character vector with features names
#' 
#' @export


select_sign_features <- function(mRList=NULL, test_method="anova", test_value = "q", cutoff_value = 0.05, feature_id="Feature_ID", values=FALSE, split_by_semicolon = FALSE) {
  
  #data.use <- match.arg(data.use, c("data", "data_ann"))
  #data <- mRList[[data.use]]
  
  stopifnot(feature_id %in% colnames(mRList$data_ann))  
  
  if(values){
    
    ans <- setNames(mRList[[test_method]][, test_value], rownames(mRList[[test_method]]))
    ans <- sort(ans[ans < cutoff_value])
    
    if(feature_id == "Name"){
      
      ans <- merge(data.frame(Feature_ID=names(ans), as.numeric(ans), stringsAsFactors = F), mRList$data_ann[, c("Feature_ID", feature_id)], by="Feature_ID")
      colnames(ans)[2] <- test_value
      
      if (split_by_semicolon) {
        if(any(grepl(";", ans[, feature_id]))){
          unl <- strsplit(ans[, feature_id], ";")
          unl <- data.frame(old=rep(ans[, feature_id], times=lengths(unl)), unlist(unl), stringsAsFactors = F)
          ans <- merge(ans, unl, by.x=3, by.y=1)
          ans <- ans[!is.na(ans[, 4]), ]
          ans <- ans[order(ans[, colnames(ans) == test_value]), ]
          idx_keep <- !duplicated(ans[, 4])
          ans <- setNames(ans[idx_keep, test_value], ans[idx_keep, 4])
        }
      }
    } else if (feature_id != "Name" & feature_id != "Feature_ID") {
      stop('feature_id must be "Feature_ID" or "Name"')
    }
    
    
  }else{
    
    mRList[[test_method]] <- mRList[[test_method]][order(mRList[[test_method]][ , test_value]),]
    ans <- row.names(mRList[[test_method]])[which(mRList[[test_method]][ , test_value] < cutoff_value)]
    
    
    if(feature_id == "Name"){
      
      ans <- as.character(mRList$data_ann[mRList$data_ann$Feature_ID %in% ans, feature_id])
      
      # 1 -> n mapping results in multiple identifiers for each feature
      
      if (split_by_semicolon) {
        if(any(grepl(";", ans))){
          ans <- unlist(strsplit(ans, ";"))
        }
      }
      
      ans <- unique(ans[!is.na(ans)])
      
    } else if (feature_id != "Name" & feature_id != "Feature_ID") {
      stop('feature_id must be "Feature_ID" or "Name"')
    }
    
  }
  
  return(ans)
  
}

