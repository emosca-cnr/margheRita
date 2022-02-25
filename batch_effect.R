#' Remove batch effect
#' @import limma
#' @export


batch_effect <- function(m_list){

  m_list$data <- limma::removeBatchEffect(m_list$data, m_list$sample_ann$batch)

  return(m_list)

}



