#' m_z_filtering
#' @param m_list margheRita mRList
#' @export

m_z_filtering <- function(m_list, do_plot=TRUE, lower_quality_mass_acc=0.4, upper_quality_mass_acc=0.8, color="black" ,...){

  m_list$metab_ann$quality <- m_list$metab_ann$mz %% 1
  all <- length(m_list$metab_ann$quality)
  idx_keep <- m_list$metab_ann$quality < lower_quality_mass_acc | m_list$metab_ann$quality > upper_quality_mass_acc

  cat("Metabolites with appropriate m/z values\n")
  cat(as.data.frame(table(idx_keep))[1, 2])
  cat("\n")
  cat("Metabolites witout appropriate m/z values\n")
  cat(all-(as.data.frame(table(idx_keep)))[1, 2])

  m_list$data <- m_list$data[idx_keep, ]
  m_list$metab_ann <- m_list$metab_ann[idx_keep, ]

  m_list$QC <- m_list$QC[idx_keep, ]

   if(do_plot){
    plot(density(m_list$metab_ann$quality), xlab="m/z Quality", main="Distrubution of m/z Quality", ...)+
    polygon(density(m_list$metab_ann$quality), col=color, border=color)
  }
  return(m_list)
}
