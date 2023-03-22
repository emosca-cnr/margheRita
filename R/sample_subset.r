#' sample_subset
#' @param mRList mRList object
#' @param samples_to_remove character vector of samples to be removed
#' @param  make_new_object if set true a new mRList will be returned
#' @return filtered mRList
#' @export
#' @return mRList_subset mRList object subset
#' @return mRList mRList object

sample_subset <- function(mRList, samples_to_remove){
   samples_to_keep <- setdiff(names(mRList$data), samples_to_remove)
    mRList$data <- mRList$data[samples_to_keep]
    mRList$data <- mRList$data[rowSums(mRList$data) > 0,]
    met_index <- as.character(rownames(mRList$data))
    mRList$metab_ann <- mRList$metab_ann[met_index,]
    mRList$sample_ann <- mRList$sample_ann[mRList$sample_ann$id %in% samples_to_keep]
    mRList$QC <- mRList$QC[met_index,]
    return(mRList)
 }