#' Remove samples from mRList
#' @param mRList object
#' @param ids rownames of $sample_ann to be removed. Otherwise, if "column" is given, all samples corresponding to the values appearing in the specfied column of $sample_ann will be removed.
#' @param column a column of $sample_ann to be used for selecting samples to be removed.
#' @export
#' @examples
#' ##To remove all samples "Blank" of the column "class"
#' mRList <- remove_samples(mRList, ids="Blank", column = "class")
#' 
#' ##To remove two samples by id
#' mRList <- remove_samples(mRList, ids=c("x186_UR_norm_AA_mealB_t00_rtech2_RP_short_SWATH_pos", "x189_UR_norm_MM_mealC_t03_rtech2_RP_short_SWATH_pos"))

remove_samples <- function(mRList=NULL, ids=NULL, column=NULL){
	
	if(is.null(column)){
		
		idx_rows <- which(rownames(mRList$sample_ann) %in% ids)
		idx_QC_rows <- which(rownames(mRList$QC_ann) %in% ids)

	}else{
		
		idx_rows <- which(mRList$sample_ann[, column] %in% ids)
		idx_QC_rows <- which(mRList$QC_ann[, column] %in% ids)
	
	}
	
	#removing
	if(length(idx_rows)>0){
		cat("removing", rownames(mRList$sample_ann)[idx_rows], "\n")
		mRList$sample_ann <- mRList$sample_ann[-idx_rows, ]
	}
	
	if(length(idx_QC_rows)>0){
		cat("removing", rownames(mRList$QC_ann)[idx_QC_rows], "\n")
		mRList$QC_ann <- mRList$QC_ann[-idx_QC_rows, ]
	}
	
	#keep data columns for the samples that are in the corresponding annotation
	mRList$data <- mRList$data[, which(colnames(mRList$data) %in% rownames(mRList$sample_ann))]
	mRList$QC <- mRList$QC[, which(colnames(mRList$QC) %in% rownames(mRList$QC_ann))]
	
	return(mRList)
	
}