#' Read files (input and metadata) both excel and .csv
#' import excel files input and metadata
#' @importFrom readxl read_excel
#' @importFrom utils read.table
#' @export

read_input_file <- function(input, metadata, split_QC=TRUE, mz_col=3, rt_col=2, data_start_col=4, MS_MS_column=NULL, type=c("excel", "txt"), ...){
  
  type <- match.arg(type)
  
  if (type=="excel"){
    data <- data.frame(readxl::read_excel(input, col_names = T), stringsAsFactors = F)
  }
  
  if (type=="txt"){
    data<-read.table(input, header=T, stringsAsFactors = F, ...)
  }
  
  mRList <- list(data=data[, -c(1:(data_start_col-1))])
  rownames(mRList$data) <- data[, 1]
  
  mRList$metab_ann <- data[, 1:(data_start_col-1)]
  colnames(mRList$metab_ann)[1] <- "Feature_ID"
  colnames(mRList$metab_ann)[mz_col] <- "mz"
  colnames(mRList$metab_ann)[rt_col] <- "rt"
  rownames(mRList$metab_ann) <- data[, 1]
  if(!is.null(MS_MS_column)){
    colnames(mRList$metab_ann)[MS_MS_column] <- "MS_MS_spectrum"
  }
  
  if (type=="excel"){
    mRList$sample_ann <- data.frame(readxl::read_excel(metadata), stringsAsFactors = F)
    rownames(mRList$sample_ann) <- mRList$sample_ann[, 1]
  }
  if (type=="txt"){
    mRList$sample_ann <- read.table(metadata, header=T, ..., stringsAsFactors = F)
    rownames(mRList$sample_ann) <- mRList$sample_ann[, 1]
  }

  if(!all(c("id", "injection_order", "batch", "class", "biological_rep", "technical_rep") %in% colnames(mRList$sample_ann))){
    cat("sample annotation must contain 'id', 'injection_order', 'batch', 'class', 'biological_rep' and 'technical_rep'\n")
    cat("found", colnames(mRList$sample_ann), "\n")
    stop("ERROR: missing mandatory sample annotation columns\n")
  }
  
  if(nrow(mRList$sample_ann) != ncol(mRList$data)){
    cat("sample annotation ", nrow(mRList$sample_ann), "\n")
    cat("data ", ncol(mRList$data), "\n")
    cat(rownames(mRList$sample_ann)[!rownames(mRList$sample_ann) %in% colnames(mRList$data)])
    cat(colnames(mRList$data)[!colnames(mRList$data) %in% rownames(mRList$sample_ann)])
    stop("ERROR: different number of elements between data and annotation")
  }
  
  idx <- match(colnames(mRList$data), mRList$sample_ann$id)
  if(any(is.na(idx))){
    print(colnames(mRList$data)[is.na(idx)])
    stop("ERROR: not all samples found in annotation")
  }
  mRList$sample_ann <- mRList$sample_ann[idx,] #ettore
  
  #SPLIT QC
  if(split_QC){
    mRList <- splitQC(mRList)
  }
  
  return(mRList)
}



