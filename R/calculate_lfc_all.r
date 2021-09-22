#' calculate_lfc_all
#' @importFrom utils combn
#' @export


calculate_lfc_all <- function(m_list, lfc_theshold=0.25) {
    samples <- unique(m_list$sample_ann$class)
    contrast <- combn(samples,2, simplify = T)
    for(i in 1:ncol(contrast)){
        sample_c <- contrast[,i]
        m_list <- calculate_lfc(m_list = m_list, contrast_samples = sample_c, lfc_theshold=0.25)
        }
    return(m_list)
}

