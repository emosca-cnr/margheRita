#' m_z_filtering
#' @import dplyr
#' @export

m_z_filtering <- function(m_list){

  #m_list$quality <- numeric(nrow(m_list$df))
  #for(i in 1:nrow(m_list$df)) {
  #  m_list$quality[i] <- (m_list$df$Average_mz[i] %% 1)
  #}

  m_list$metab_ann$quality <- m_list$metab_ann$Average_mz %% 1
  plot(density(m_list$metab_ann$quality))

  idx_keep <- m_list$metab_ann$quality < 0.4 | m_list$metab_ann$quality > 0.8

  m_list$df <- m_list$df[idx_keep, ]
  m_list$metab_ann <- m_list$metab_ann[idx_keep, ]

  #m_list$df <- dplyr::filter(m_list$df, m_list$quality < 0.4 | m_list$quality > 0.8)

  return(m_list)

}
