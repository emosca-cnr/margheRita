#' heatscatter_chromatography
#' @importFrom LSD heatscatter
#' @importFrom grDevices png
#' @export

heatscatter_chromatography <- function(mRList, mz_limits=NULL, rt_limits=NULL, sample=NULL, outfile="heatscatter.png" , ...) {
   if (!is.null(sample)) {
      df <- mRList$metab_ann[, c("rt", "mz", paste(sample[1], "_mean", sep = ""))]
      df <- df[df[,3]> 0,]
      df <- df[1:2]
   }

   df <- mRList$metab_ann[, c("rt", "mz")]
   df[is.na(df)] <- 0


   if (!is.null(mz_limits) | !is.null(rt_limits)) {
      df <- df[df$rt >= rt_limits[1] &  df$rt <= rt_limits[2], ]
      df <- df[df$mz >= mz_limits[1] &  df$mz <= mz_limits[2], ]
   }

   grDevices::png(outfile, width = 180, height = 180, res=300, units="mm")
   LSD::heatscatter(df$rt, df$mz, xlab="RT", ylab="m/z", ...)
   dev.off()

}
