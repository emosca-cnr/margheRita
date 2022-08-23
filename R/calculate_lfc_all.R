#' Calculate the log fold change between all pairs of sample classes
#' @importFrom utils combn
#' @param  mRList mRList object
#' @export
#' @return mRList object with fold changes in mRList$metab_ann 


calculate_lfc_all <- function(mRList) {
    samples <- unique(mRList$sample_ann$class)
    contrast <- combn(samples, 2, simplify = T)
    for(i in 1:ncol(contrast)){
        sample_c <- contrast[, i]
        mRList <- calculate_lfc(mRList = mRList, contrast_samples = sample_c)
    }
    return(mRList)
}

