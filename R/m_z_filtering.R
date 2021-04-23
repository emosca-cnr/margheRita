#' @import dplyr
#' @export

m_z_filtering <- function(df, m_z_average){
  df$quality <- NA
  for(i in 1:nrow(df)) {
    df$quality[i] <- (df$m_z_average[i] %% 1)
  }
  df <- filter(df, df$quality < 0.4 | df > 0.8)
  qual_plot <- plot(density(df$quality))
  return(as.data.frame(df), qual_plot)
}
