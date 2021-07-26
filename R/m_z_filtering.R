#' m_z_filtering
#' @export

m_z_filtering <- function(m_list, do_plot=TRUE){

  #m_list$quality <- numeric(nrow(m_list$df))
  #for(i in 1:nrow(m_list$df)) {
  #  m_list$quality[i] <- (m_list$df$Average_mz[i] %% 1)
  #}

  m_list$metab_ann$quality <- m_list$metab_ann$mz %% 1
  if(do_plot){
    plot(density(m_list$metab_ann$quality))
  }

  idx_keep <- m_list$metab_ann$quality < 0.4 | m_list$metab_ann$quality > 0.8

  cat("Metabolites with appropriate m/z values\n")
  print(table(idx_keep))

  m_list$data <- m_list$data[idx_keep, ]
  m_list$metab_ann <- m_list$metab_ann[idx_keep, ]

  m_list$QC <- m_list$QC[idx_keep, ]


  return(m_list)

}
