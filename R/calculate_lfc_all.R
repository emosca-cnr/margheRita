#' calculate_lfc_all
#' @importFrom utils combn
#' @export


calculate_lfc_all <- function(mRList, lfc_theshold=0.25) {
    samples <- unique(mRList$sample_ann$class)
    contrast <- combn(samples,2, simplify = T)
    for(i in 1:ncol(contrast)){
        sample_c <- contrast[,i]
        mRList <- calculate_lfc(mRList = mRList, contrast_samples = sample_c, lfc_theshold=0.25)
        }
    return(mRList)
}

