#' heatscatter_chromatography
#' @importFrom LSD heatscatter
#' @importFrom grDevices dev.off pdf
#' @export

heatscatter_chromatography <- function(m_list, mz_limits=NULL, rt_limits=NULL, ...) {

   df <- m_list$metab_ann[, c("rt", "mz")]
   df[is.na(df)] <- 0

   if (!is.null(mz_limits) | !is.null(rt_limits)) {
      df <- df[df$rt >= rt_limits[1] &  df$rt <= rt_limits[2], ]
      df <- df[df$mz >= mz_limits[1] &  df$mz <= mz_limits[2], ]
   }

   LSD::heatscatter(df$rt, df$mz, xlab="Retenction Time", ylab="m/z feature", ...)

}
