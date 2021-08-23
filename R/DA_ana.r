#' DA_ana
#' @export

DA_ana <- function(m_list, FC_theshold, contrast_samples){
    log_data <- m_list$metab_ann[grep(paste(contrast_samples[1], "_mean", "|",contrast_samples[2], "_mean", sep = ""), names(m_list$metab_ann))]
    log_data[log_data == 0 | is.na(log_data), ] <- 1
    log_data <- log2(log_data)
    log_data$Log2FC <- log_data[2] - log_data[1]
    new_name <- paste("Log2_FC",contrast_samples[1], contrast_samples[2], sep="_")
    names(log_data)[3] <- new_name
    m_list$metab_ann <- cbind(m_list$metab_ann, log_data[3])
    return(m_list)
}