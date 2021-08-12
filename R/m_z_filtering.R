#' m_z_filtering
#' @export

m_z_filtering <- function(m_list, do_plot=TRUE, lower_quality_mass_acc, upper_quality_mass_acc){

  m_list$metab_ann$quality <- m_list$metab_ann$mz %% 1

  idx_keep <- m_list$metab_ann$quality < lower_quality_mass_acc | m_list$metab_ann$quality > upper_quality_mass_acc

  cat("Metabolites with appropriate m/z values\n")
  print(table(idx_keep))

  m_list$data <- m_list$data[idx_keep, ]
  m_list$metab_ann <- m_list$metab_ann[idx_keep, ]

  m_list$QC <- m_list$QC[idx_keep, ]

   if(do_plot){
    plot(density(m_list$metab_ann$quality))
  }

  return(m_list)

}
