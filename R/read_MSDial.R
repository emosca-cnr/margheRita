#' Read MSDial file header and prepare a sample annotation file
#' @description Import data and metadata as excel or .csv or .txt format
#' @param MSDial_file MSDial file name
#' @param nrows to extract MSDial sample info only the top nrows will be read
#' @param ... further arguments to read.delim
#' @importFrom utils read.delim write.table
#' @return data.frame with sample info
#' @export

read_MSDial <- function(MSDial_file=NULL, nrows=50, ...){
  
  
  cat("Reading MS-Dial input file...\n")
  feat_data <- read.delim(MSDial_file, stringsAsFactors = F, quote="", comment.char="", header=F, nrows = nrows, ...)
  
  feat_data_header_row <- which(grepl("Alignment ID", feat_data[, 1]))
  sample_data_header_col_start <- which(grepl("Class", feat_data[1, ]))
  sample_data_header_row_end <-  which(feat_data[, sample_data_header_col_start] == "MS/MS spectrum")
  
  sample_col <- which(!is.na(feat_data[2, ]))[-c(1:sample_data_header_col_start)] #indices of columns that contain samples
  data_end_col <- sample_col[length(sample_col)]
  
  
  sample_ann <- t(feat_data[1:(sample_data_header_row_end), sample_data_header_col_start:data_end_col])
  colnames(sample_ann) <- sample_ann[1, ]
  sample_ann <- data.frame(sample_ann[-1, ])
  sample_ann <- sample_ann[sample_ann$MS.MS.spectrum != 1, ] #remove average and std
  colnames(sample_ann)[colnames(sample_ann) == "MS.MS.spectrum"] <- "id"
  colnames(sample_ann)[colnames(sample_ann) == "Injection.order"] <- "injection_order"
  colnames(sample_ann)[colnames(sample_ann) == "Batch.ID"] <- "batch"
  colnames(sample_ann)[colnames(sample_ann) == "Class"] <- "class"
  sample_ann$biological_rep <- sample_ann$technical_rep <- ""
  
  sample_ann <- sample_ann[, c("id", "injection_order", "batch", "class", "biological_rep", "technical_rep")]
  
  return(sample_ann)
  
}



