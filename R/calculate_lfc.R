#' calculate_lfc
#' @export

calculate_lfc <- function(m_list, lfc_theshold=0.25, contrast_samples){

    log_data <- m_list$metab_ann[, grep(paste(contrast_samples[1], "_mean", "|",contrast_samples[2], "_mean", sep = ""), names(m_list$metab_ann))]

    log_data[is.na(log_data)] <- 0
    log_data <- log2(log_data+1)
    log_data$Log2FC <- log_data[, 2] - log_data[, 1]

    names(log_data)[3] <- paste("Log2_FC",contrast_samples[1], contrast_samples[2], sep="__")

    m_list$metab_ann <- cbind(m_list$metab_ann, log_data[, 3])
    names(m_list$metab_ann)[length(names(m_list$metab_ann))] <- names(log_data)[3]
    m_list$metab_ann[length(names(m_list$metab_ann))]
    m_list$metab_ann[,ncol(m_list$metab_ann)][m_list$metab_ann[,ncol(m_list$metab_ann)] > -(lfc_theshold) & m_list$metab_ann[,ncol(m_list$metab_ann)] < lfc_theshold] <- 0


    return(m_list)
}
