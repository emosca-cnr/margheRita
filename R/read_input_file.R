#' Create an mRList from feature data and sample information
#' @description Import data and metadata as excel or .csv or .txt format
#' @param split_QC (default=TRUE) split input and metadata in two one for QC and the other for data
#' @param MSDial (default=TRUE) the input file is MS-Dial native format
#' @param mz_col column that contains m/z value
#' @param rt_col column that contains rt value
#' @param data_start_col column from which data starts
#' @importFrom utils read.delim
#' @export
#' @param MS_MS_column column that contains the MS/MS spectrum
#' @param ... further arguments to read.delim
#' @param feature_file tab-delimited text files with feature abundances, mz, rt and MS/MS spectra
#' @param sample_file sample information file
#' @return mRList object

read_input_file <- function(feature_file=NULL, sample_file=NULL, MSDial=TRUE, split_QC=TRUE, mz_col=NULL, rt_col=NULL, data_start_col=NULL, MS_MS_column=NULL, ...){


  if(MSDial){

    cat("Reading MS-Dial input file...\n")
    feat_data <- read.delim(feature_file, stringsAsFactors = F, quote="", comment.char="", header=F)

    stopifnot(any(feat_data[2, ] == "File type"), any(feat_data[4, ] == "Average"))
    data_start_col <- which(feat_data[2, ] == "File type")+1
    data_end_col <- which(feat_data[4, ] == "Average")[1]-1
    if(is.na(data_end_col)){
      data_end_col <- ncol(feat_data)
    }
    
    stopifnot(any(feat_data[5, ] == "Average Mz"), any(feat_data[5, ] == "Average Rt(min)"))
    mz_col <- which(feat_data[5, ] == "Average Mz")
    rt_col <- which(feat_data[5, ] == "Average Rt(min)")
    MS_MS_column <- which(feat_data[5, ] == "MS/MS spectrum")
    Ontology_column <- which(feat_data[5, ] == "Ontology")

    feat_data <- feat_data[-c(1:4), ]
    colnames(feat_data) <- feat_data[1, ]
    feat_data <- feat_data[-1, ]
    rownames(feat_data) <- feat_data[, 1] <- paste0("F", feat_data[, 1])

    stopifnot(any(colnames(feat_data) == "Metabolite name"), any(colnames(feat_data)== "SMILES"))
    mRList <- list(
      data=data.frame(feat_data[, data_start_col:data_end_col]),
      metab_ann=data.frame(feat_data[, c(1, which(colnames(feat_data) == "Metabolite name"), which(colnames(feat_data) == "SMILES"), Ontology_column, rt_col, mz_col, MS_MS_column)])
    )
    colnames(mRList$metab_ann) <- c("Feature_ID", "MSDialName", "MSDialSMILES", "Ontology", "rt", "mz", "MS_MS_spectrum")
    mRList$metab_ann$rt <- as.numeric(mRList$metab_ann$rt)
    mRList$metab_ann$mz <- as.numeric(mRList$metab_ann$mz)

    for(i in 1:ncol(mRList$data)){
      mRList$data[, i] <- as.numeric(mRList$data[, i])
    }

  }else{

    cat("Reading generic input file...\n")

    data <- read.delim(feature_file, header=T, stringsAsFactors = F, quote="", comment.char="")

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
  }

  print(dim(mRList$data))

  cat("Reading sample annotation\n")
  mRList$sample_ann <- read.delim(sample_file, header=T, stringsAsFactors = F, ...)
  mRList$sample_ann$id <- make.names(mRList$sample_ann$id)
  #rownames(mRList$sample_ann) <- mRList$sample_ann[, 1]
  rownames(mRList$sample_ann) <- mRList$sample_ann$id
  
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

  # idx <- match(colnames(mRList$data), mRList$sample_ann$id)
  # 
  # if(any(is.na(idx))){ #GF: fixing special characters
  #   colnames(mRList$data) <- gsub("[[:punct:]]", ".", colnames(mRList$data))
  #   colnames(mRList$data) <- gsub("\\s", ".", colnames(mRList$data))
  #   colnames(mRList$data) <- ifelse(grepl("^\\d", colnames(mRList$data)), paste0("X", colnames(mRList$data)), colnames(mRList$data))
  # 
  #   mRList$sample_ann$id <- gsub("[[:punct:]]", ".", mRList$sample_ann$id)
  #   mRList$sample_ann$id <- gsub("\\s", ".", mRList$sample_ann$id)
  #   mRList$sample_ann$id <- ifelse(grepl("^\\d", mRList$sample_ann$id), paste0("X", mRList$sample_ann$id), mRList$sample_ann$id)
  # }
  
  idx_updt <- match(colnames(mRList$data), mRList$sample_ann$id)

  if(any(is.na(idx_updt))){
    print(colnames(mRList$data)[is.na(idx_updt)])
    stop("ERROR: not all samples found in annotation")
  }
  
  mRList$sample_ann <- mRList$sample_ann[idx_updt,] #ettore


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



