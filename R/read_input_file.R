#' Import data and metadata
#' @description Import data and metadata as excel or .csv or .txt format
#' @param split_QC (default=TRUE) split input and metadata in two one for QC and the other for data
#' @param  mz_col column that contains m/z value
#' @param rt_col column that contains rt value
#' @param data_start_col column from which data starts
#' @param type type of imported files Excel or .csv or .txt
#' @importFrom readxl read_excel
#' @importFrom utils read.table
#' @export
#' @param MS_MS_column column that contains the MS/MS spectrum
#' @param ... further arguments to read.table
#' @param input data file name
#' @param metadata metadata file name
#' @return mRList S4 object that contain data and metadata.
#' @examples
#' ##library(dataset.margheRita)
#' ##dataset(norm_pos)
#' mRList<-(norm_pos.xlsx, norm_pos_meta.xlsx, split_QC=TRUE, mz_col=3, rt_col=2, data_start_col=4, MS_MS_column=NULL, type="Excel")

read_input_file <- function(input=NULL, metadata=NULL, split_QC=TRUE, mz_col=3, rt_col=2, data_start_col=5, MS_MS_column=NULL, type=c("Excel", "txt"), ...){

  type <- match.arg(type)

  if (type=="Excel"){
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

  if (type=="Excel"){
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


  ###
  cat("# Samples:", nrow(mRList$sample_ann), "\n")

  #SPLIT QC
  if(split_QC){
    mRList <- splitQC(mRList)
    cat("# QC Samples:", nrow(mRList$QC_ann), "\n")
  }

  #class by biological replicates
  cat("Class by biological replicates\n")
  print(table(mRList$sample_ann$class, mRList$sample_ann$biological_rep))


  return(mRList)
}



