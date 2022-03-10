#' calculate_lfc
#' @export

calculate_lfc <- function(mRList, lfc_theshold=0.25, contrast_samples){

    log_data <- mRList$metab_ann[, grep(paste(contrast_samples[1], "_mean", "|",contrast_samples[2], "_mean", sep = ""), names(mRList$metab_ann))]

    log_data[is.na(log_data)] <- 0
    log_data <- log2(log_data+1)
    log_data$Log2FC <- log_data[, 2] - log_data[, 1]

    names(log_data)[3] <- paste("Log2_FC",contrast_samples[1], contrast_samples[2], sep="__")

    mRList$metab_ann <- cbind(mRList$metab_ann, log_data[, 3])
    names(mRList$metab_ann)[length(names(mRList$metab_ann))] <- names(log_data)[3]
    mRList$metab_ann[length(names(mRList$metab_ann))]
    mRList$metab_ann[,ncol(mRList$metab_ann)][mRList$metab_ann[,ncol(mRList$metab_ann)] > -(lfc_theshold) & mRList$metab_ann[,ncol(mRList$metab_ann)] < lfc_theshold] <- 0


    return(mRList)
}
