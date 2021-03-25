#' pareto scaling of the dataframe after missing value correction
# dataframe only with samples and intensity of each m/z and RT
#' @param df data.frame without the first three columns (MS DIAL, RT and m/z)
#' @importFrom utils write.csv
#' @export
#' @importFrom stats sd


pareto <- function(df) {
    # Here we extract numeric data and perform Pareto scaling
    sample_classes <- df[, 1:3]
    x <- df[, 4:dim(df)[2]]
   # Here we perform centering
  x.centered <- apply(x, 1, function(x) x - mean(x))
  # Then we perform scaling on the mean-centered matrix
  x.sc <- apply(x.centered, 1, function(x) x/sqrt(sd(x)))
  x.sc <- cbind(sample_classes, x.sc)
  utils::write.csv(x.sc, "paretoscaled.csv",sep=",")
  return(x.sc)
}

