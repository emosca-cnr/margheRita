#' Load nutrition dataset
#' @param column POS or NEG
#' @param normalized dataset with/without normalized samples.#' 
#' @export

load_nutrition_dataset <- function(column=NULL, normalized=FALSE){
	
	column <- match.arg(column, c("POS", "NEG"))
	
	cat("Column", column, "\n")
	cat("Normalized", normalized, "\n")
	
	input_data_file <- paste0("nutrition_", ifelse(normalized, "norm_", ""), column, ".rds")
	input_data_file <- system.file("extdata", input_data_file, package = "margheRita")
	
	mRList <- readRDS(input_data_file)
	
	return(mRList)

}