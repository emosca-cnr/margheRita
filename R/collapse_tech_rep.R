#' Collapse technical replicates
#'
#' @description Collapse technical replicates of each biological replicate calculating the mean.
#' The resulting dataframe and metadata contain only biological replicates.
#'
#' @param mRList mRList object
#' @param remove.QC whether to remove QC samples
#' @export
#' @return mRList object with only biological replicates

collapse_tech_rep <-function(mRList=NULL, remove.QC=TRUE){

  cat("According to dataset size, this might take a few minutes.\n")
  mRList$sample_ann$class_biorep <- as.factor(paste(mRList$sample_ann$class, mRList$sample_ann$biological_rep, sep="_"))

  mRList$data <- t(apply(mRList$data, 1, function(x) tapply(x, mRList$sample_ann$class_biorep, mean)))
  
  mRList$sample_ann <- unique(mRList$sample_ann[, c("class_biorep", "class", "biological_rep")])
  rownames(mRList$sample_ann) <- mRList$sample_ann$class_biorep


  if(remove.QC){
    mRList$QC <- mRList$QC_ann <- NULL
  }else{
    mRList$QC_ann$class_biorep <- as.factor(paste(mRList$QC_ann$class, mRList$QC_ann$biological_rep, sep="_"))
    mRList$QC <- t(apply(mRList$QC, 1, function(x) tapply(x, mRList$QC_ann$class_biorep, mean)))
    mRList$QC_ann <-  unique(mRList$QC_ann[, c("class_biorep", "class", "biological_rep")])
    rownames(mRList$QC_ann) <- mRList$QC_ann$class_biorep
    
  }

  #ensure correct order of samples
  idx <- match(colnames(mRList$data), mRList$sample_ann$class_biorep)
  if(any(is.na(idx))){
    stop("ERROR: not all samples found in annotation")
  }
  mRList$sample_ann <- mRList$sample_ann[idx,] #ettore

  if(nrow(mRList$sample_ann) != ncol(mRList$data)){
    cat("data ", nrow(mRList$sample_ann), "\n")
    cat("data ", ncol(mRList$data), "\n")
    stop("ERROR: different number of elements between data and annotation")
  }

  #ensure correct order of metabolites
  if(!identical(rownames(mRList$data), mRList$metab_ann[, 1])){
    idx <- match(rownames(mRList$data), mRList$metab_ann[, 1])
    if(any(is.na(idx))){
      stop("ERROR: not all metabolites found.")
    }
    mRList$metab_ann <- mRList$metab_ann[idx, ] #ettore
  }

  return(mRList)
}
