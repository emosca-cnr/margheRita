#' Keep the feature with the highest average across samples
#' @param X features-by-samples data matrix
#' @export
#' 
keep_strongest_representative <- function(X=NULL){
  
  #calculate rowMeans
  rn_avg <- cbind(X, avg=rowMeans(X[, -c(1:2)]))
  
  #sort by decreasing rowMeans
  rn_avg <- rn_avg[order(-rn_avg$avg), ]
  
  #find duplicates among novel row.names
  idx_dupl <- which(duplicated(rn_avg$Name))
  
  if(length(idx_dupl)>0){
    
    rn_avg <- rn_avg[-idx_dupl, ]
    
    #if there are NA, remove the row
    idx_NA <- which(is.na(rn_avg$Name))
    if(length(idx_NA)>0){
      rn_avg <- rn_avg[-idx_NA, ]
    }
    
    rownames(rn_avg) <- rn_avg$Name
    rn_avg$avg <- NULL
    rn_avg <- rn_avg[order(rownames(rn_avg)), ]
    
  }else{
    
    cat("No duplicated features.\n")
    rn_avg <- NULL
    
  }
  
  return(rn_avg)
  
}