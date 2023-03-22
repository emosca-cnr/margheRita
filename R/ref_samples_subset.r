#' remove_ref_samples
#' @param mRList mRList object
#' @param samples_to_remove character vector of samples to be removed
#' @param  make_new_object if set true a new mRList will be returned
#' @param QCs_to_remove vector of QCs names samples to remove
#' @return filtered mRList
#' @export
#' @return mRList_subset mRList object subset
#' @return mRList mRList object
#' @import data.table %in%

remove_ref_samples <- function(mRList, remove_all_QCs=FALSE, remove_all_Blanch=FALSE, QCs_to_remove=NULL){
    if(remove_all_Blanch){
        samples_to_keep <- mRList$sample_ann[mRList$sample_ann$class != "Blank",][,1]
        mRList$data <- mRList$data[samples_to_keep]
        mRList$data <- mRList$data[rowSums(mRList$data) > 0,]
        met_index <- as.character(rownames(mRList$data))
        mRList$metab_ann <- mRList$metab_ann[met_index,]
        mRList$sample_ann <- mRList$sample_ann[mRList$sample_ann$id %in% samples_to_keep,]
    }
    if(remove_all_QCs){
        mRList <- mRList[1:3]
    }
    if(!is.null(QCs_to_remove)){
        samples_to_keep <- setdiff(names(mRList$QC), QCs_to_remove)
        mRList$QC <- mRList$QC[samples_to_keep]
        mRList$QC <- mRList$QC[rowSums(mRList$QC) > 0,]
        met_index <- as.character(rownames(mRList$QC))
        mRList$QC <- mRList$QC[met_index,]
        mRList$data <- mRList$data[met_index,]
        mRList$metab_ann <- mRList$metab_ann[met_index,]
        mRList$QC_ann <- mRList$QC_ann[mRList$QC_ann$id %in% samples_to_keep]
    }
    return(mRList)
 }