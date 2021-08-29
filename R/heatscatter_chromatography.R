#' heatscatter_chromatography
#' @importFrom LSD heatscatter
#' @importFrom grDevices dev.off pdf
#' @export

heatscatter_chromatography <- function(m_list, mz_limits=NULL, rt_limits=NULL, stat=F, color_palette=1, title=NULL) {
   df <- as.data.frame(m_list$metab_ann)[2:3]
   df[is.na(df),] <- 0
   if (mz_limits!=NULL | rt_limits!=NULL) {
      df <- df[df$rt >= rt_limits[1] &  df$rt <= rt_limits[2],]
      df <- df[df$mz >= mz_limits[1] &  df$mz <= mz_limits[2],]
   }
   heatscatter(df$rt, df$mz, 
      xlab="Retenction Time",
      ylab="m/z feature",
      simulate=ifelse(color_palette==1, F, T),
      main=ifelse(title==NULL, NULL, title),
      add.contour=ifelse(stat==F, F, T))
      return(m_list)
}