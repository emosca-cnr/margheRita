#' Draw the "heatscatter chromatography"
#' @importFrom LSD heatscatter
#' @importFrom grDevices png
#' @export
#' @param mRList mRList object
#' @param mz_limits numeric vector with minimum and maximum rt values
#' @param rt_limits numeric vector with minimum and maximum rt values
#' @param sample ???
#' @param outfile out file name
#' @param ... further arguments to LSD::heatscatter
#'

heatscatter_chromatography <- function(mRList, mz_limits=NULL, rt_limits=NULL, sample=NULL, outfile="heatscatter.png" , ...) {

  ### change sample into colnames or index
  if (!is.null(sample)) {
    col_index <- grep(sample, names(mRList$metab_ann))
    col_names <- append(c("rt", "mz"), names(mRList$metab_ann)[col_index])
    df <- mRList$metab_ann[col_names]
    df$sum <- rowSums(df[3:ncol(df)])
    df <- df[df[,ncol(df)]> 0,]
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
