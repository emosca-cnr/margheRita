#' heatscatter_chromatography
#' @importFrom LSD heatscatter
#' @importFrom grDevices dev.off pdf
#' @export

heatscatter_chromatography <- function(m_list, mz_limits, rt_limits, stat, color_palette, title) {
   if (mz_limits==NULL | rt_limits==NULL) {
      rt <- m_list$metab_ann$rt
      mz <- m_list$metab_ann$mz
   } else {
       rt <- m_list$metab_ann$rt[m_list$metab_ann$rt > rt_limits[1]  & m_list$metab_ann$rt < rt_limits[2]]
       mz <- m_list$metab_ann$mz[m_list$metab_ann$mz > mz_limits[1]  & m_list$metab_ann$mz < mz_limits[2]]
   }
   heatscatter(m_list$metab_ann$rt, m_list$metab_ann$mz, 
      xlab="Retenction Time",
      ylab="m/z feature",
      simulate=ifelse(color_palette==1, F, T),
      main=ifelse(title==NULL, NULL, title),
      add.contour=ifelse(stat==F, F, T))
return(m_list)
}