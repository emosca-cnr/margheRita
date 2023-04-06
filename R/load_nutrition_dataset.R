#' Load nutrition dataset
#' @param column POS or NEG
#' @param normalized dataset with/without experimental normalization to Urine volume 
#' @export

load_nutrition_dataset <- function(mode=NULL, normalized=FALSE){
	
	mode <- match.arg(mode, c("POS", "NEG"))
	
	cat("mode", mode, "\n")
	cat("Normalized", normalized, "\n")
	
	input_data_file <- paste0("nutrition_", ifelse(normalized, "norm_", ""), mode, ".rds")
	input_data_file <- system.file("extdata", input_data_file, package = "margheRita")
	
	mRList <- readRDS(input_data_file)
	
	return(mRList)

}