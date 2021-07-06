#' Read files (input and metadata) both excel and .csv, default .csv

#' import excel files input and metadata
#' @import readxl
#' @export

read_input_file <- function(input, metadata, split_QC=TRUE){

  data <- data.frame(readxl::read_excel(input, col_names = T), stringsAsFactors = F)

  m_list <- list(data=data[, -c(1:3)])
  rownames(m_list$data) <- data[, 1]

  m_list$metab_ann <- data[, 1:3]
  rownames(m_list$metab_ann) <- data[, 1]

  m_list$sample_ann <- data.frame(readxl::read_excel(metadata), stringsAsFactors = F)
  rownames(m_list$sample_ann) <- m_list$sample_ann[, 1]

  if(!all(c("id", "injection_order", "batch", "class", "biological_rep", "technical_rep") %in% colnames(m_list$sample_ann))){
    cat("sample annotation must contain 'id', 'injection_order', 'batch', 'class', 'biological_rep' and 'technical_rep'\n")
    cat("found", colnames(m_list$sample_ann), "\n")
    stop("ERROR: missing mandatory sample annotation columns\n")
  }

  if(nrow(m_list$sample_ann) != ncol(m_list$data)){
    cat("data ", nrow(m_list$sample_ann), "\n")
    cat("data ", ncol(m_list$data), "\n")
    stop("ERROR: different number of elements between data and annotation")
  }

  idx <- match(colnames(m_list$data), m_list$sample_ann$id)
  if(any(is.na(idx))){
    stop("ERROR: not all samples found in annotation")
  }
  m_list$sample_ann <- m_list$sample_ann[idx,] #ettore

  #SPLIT QC
  if(split_QC){
    m_list <- splitQC(m_list)
  }

  return(m_list)
}

#import csv files MSdial ID as rownames(metabolite)

# read_input_file <- function(input,metadata){
# df <- read.delim(file,
#                  sep = ",",
#                  h = T,
#                  row.names = 1)
# metadata <- read.delim(metadata,
#                        sep = ",",
#                        h = T,
#                        row.names = 1)
# }

#setmetabolites name
#metabolites <- dataframe[, c(1:3)] #ID MSDial+RT+average_m/z

#Create an additional column called "Group" in the metadata in which
#according to subclass (diet, genotype, etc.) a number as class is associated
#Code is specific for Metadata used for example, if it is useful it can be generalized


#subclass <- metadata$subclass
#print(subclass)
#class <- metadata$class
#print(class)

#@import library(dplyr)

#group <- function(metadata) {
#metadata <- metadata %>%
# mutate(
#Group = case_when(
# metadata$subclass == "STAND" ~ "1",
# metadata$subclass == "OMEGA" ~ "2",
# metadata$class == "QC" ~ "3",
#)
# )
# metadata$Group <- as.factor(metadata$Group)
#}

