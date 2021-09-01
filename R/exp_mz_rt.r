#' exp_mz_rt
#' @importFrom plotly plot_ly
#' @export

exp_mz_rt <- function(m_list, sample=c(), mz_limits=NULL, rt_limits=NULL) {
   df <- m_list$metab_ann[,c("rt", "mz", names(m_list$metab_ann)[grep(paste(sample[1], "_mean", sep = ""), colnames(m_list$metab_ann))], "quality")]
   df[is.na(df)] <- 0
   if (!is.null(mz_limits) | !is.null(rt_limits)) {
      df <- df[df$rt >= rt_limits[1] & df$rt <= rt_limits[2],]
      df <- df[df$mz >= mz_limits[1] & df$mz <= mz_limits[2],]
   }
   plotly::plot_ly(df, x= ~rt, y=~mz, z= df[,3], type = 'scatter3d',
        marker = list(color = ~quality, colorscale = c('#FFE1A1', '#683531'),  
                      symbol = "circle", 
                      size = 2, alpha = 0.1,
                      showscale = TRUE))
}
