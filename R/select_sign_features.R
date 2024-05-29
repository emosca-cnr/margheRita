#' select significant features
#' 
#' The output of this function can be passed in the features argument of the h_map function
#' 
#' @param mRList mRList object
#' @param data.use whether to use mRList$data (default) or mRList$data_ann, which is the annotated version obtained via metabolite_identification
#' @param test_method the name of the statistical test table contained in the mRlist object
#' @param test_value the column of the statistical test table to consider (e.g.: p or q)
#' @param cutoff_value the value of p-value cut-off to consider a feature as significant
#' 
#' #' @return a character vector with features names
#' 
#' @export
select_sign_features <- function(mRList=NULL, data.use = c("data", "data_ann"), test_method="anova", test_value = "q", cutoff_value = 0.05) {
  
  data.use <- match.arg(data.use, c("data", "data_ann"))
  
  data <- mRList[[data.use]]
  
  mRList[[test_method]] <- mRList[[test_method]][order(mRList[[test_method]][ , test_value]),]
  
  sign_feat <- row.names(mRList[[test_method]])[which(mRList[[test_method]][ , test_value] < cutoff_value)]
  
  
  
  if (data.use == "data") {
    
    return(sign_feat[which(sign_feat %in% row.names(data))])
    
  } else if (data.use == "data_ann") {
    
    
    mRList[["metab_ann"]] <- mRList[["metab_ann"]][match(row.names(mRList[[test_method]]), mRList[["metab_ann"]]$Feature_ID), ]
    
    sign_molecules_all <- mRList[["metab_ann"]][which(mRList[["metab_ann"]]$Feature_ID %in% sign_feat), "Name"]
    
    sign_molecules <- unique(sign_molecules_all[which(!is.na(sign_molecules_all))])
    
    return(row.names(data)[which(row.names(data) %in% sign_molecules)])
    
  }
}

