#' Collapse technical replicates of each biological replicate of samples, calculating the mean.
#' @param mRList mRList
#' @importFrom stats aggregate
#' @export
#' @return mRList with only biological replicates


collapse_tech_rep <-function(mRList, remove.QC=FALSE){

  cat("According to dataset size, this might take a few minutes.\n")
  mRList$sample_ann$class_biorep <- as.factor(paste(mRList$sample_ann$class, mRList$sample_ann$biological_rep, sep="_"))

  mRList$data <- stats::aggregate.data.frame(t(mRList$data), list(mRList$sample_ann$class_biorep), mean)

  rownames(mRList$data) <- mRList$data$Group.1
  mRList$data <- mRList$data[, -1]
  mRList$data <- t(mRList$data)
  mRList$data <- as.data.frame(mRList$data)
  #row.names(mRList$data)<-mRList$metab_ann$MS.Dial.ID

  mRList$sample_ann <- unique(mRList$sample_ann[, c("class_biorep", "class", "biological_rep")])
  rownames(mRList$sample_ann) <- mRList$sample_ann$class_biorep


  if(remove.QC){
    mRList$QC <- mRList$QC_ann <- NULL
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
