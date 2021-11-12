#' Collapse technical replicates of each biological replicate of samples, calculating the mean.
#' @param m_list mRList
#' @importFrom stats aggregate
#' @export
#' @return m_list with only biological replicates


collapse_tech_rep <-function(m_list, remove.QC=FALSE){

  cat("According to dataset size, this might take a few minutes.\n")
  m_list$sample_ann$class_biorep <- as.factor(paste(m_list$sample_ann$class, m_list$sample_ann$biological_rep, sep="_"))

  m_list$data <- stats::aggregate.data.frame(t(m_list$data), list(m_list$sample_ann$class_biorep), mean)

  rownames(m_list$data) <- m_list$data$Group.1
  m_list$data <- m_list$data[, -1]
  m_list$data <- t(m_list$data)
  m_list$data <- as.data.frame(m_list$data)
  #row.names(m_list$data)<-m_list$metab_ann$MS.Dial.ID

  m_list$sample_ann <- unique(m_list$sample_ann[, c("class_biorep", "class", "biological_rep")])
  rownames(m_list$sample_ann) <- m_list$sample_ann$class_biorep


  if(remove.QC){
    m_list$QC <- m_list$QC_ann <- NULL
  }

  #ensure correct order of samples
  idx <- match(colnames(m_list$data), m_list$sample_ann$class_biorep)
  if(any(is.na(idx))){
    stop("ERROR: not all samples found in annotation")
  }
  m_list$sample_ann <- m_list$sample_ann[idx,] #ettore

  if(nrow(m_list$sample_ann) != ncol(m_list$data)){
    cat("data ", nrow(m_list$sample_ann), "\n")
    cat("data ", ncol(m_list$data), "\n")
    stop("ERROR: different number of elements between data and annotation")
  }

  #ensure correct order of metabolites
  if(!identical(rownames(m_list$data), m_list$metab_ann[, 1])){
    idx <- match(rownames(m_list$data), m_list$metab_ann[, 1])
    if(any(is.na(idx))){
      stop("ERROR: not all metabolites found.")
    }
    m_list$metab_ann <- m_list$metab_ann[idx, ] #ettore
  }

  return(m_list)
}
