#' exclude features with problematic m/z values
#' @import dplyr
#' @export

m_z_filtering <- function(df, m_z_average){
  df$quality <- NA
  for(i in 1:nrow(df)) {
    df$quality[i] <- (df$Average_mz[i] %% 1)
  }
  filter(df, df$quality < 0.4 | df > 0.8)
  plot(density(df$quality))
}
