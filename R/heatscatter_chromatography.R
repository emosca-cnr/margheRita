#' heatscatter_chromatography
#' @importFrom LSD heatscatter
#' @importFrom grDevices dev.off pdf
#' @export

heatscatter_chromatography <- function(m_list, mz_limits=NULL, rt_limits=NULL, sample=NULL,...) {
   if (!is.null(sample)) {
      df <- m_list$metab_ann[, c("rt", "mz", paste(sample[1], "_mean", sep = ""))]
      df <- df[df[,3]> 0,]
      df <- df[1:2]
   }

   df <- m_list$metab_ann[, c("rt", "mz")]
   df[is.na(df)] <- 0
   

   if (!is.null(mz_limits) | !is.null(rt_limits)) {
      df <- df[df$rt >= rt_limits[1] &  df$rt <= rt_limits[2], ]
      df <- df[df$mz >= mz_limits[1] &  df$mz <= mz_limits[2], ]
   }

   LSD::heatscatter(df$rt, df$mz, xlab="Retenction Time", ylab="m/z feature", ...)

}
