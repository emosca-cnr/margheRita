#' Plot the distribution of m/z decimal parts
#' @param mRList mRList object
#' @importFrom ggplot2 ggplot geom_histogram ylab xlab set_theme
#' @return ggplot object
#' @export
#' 
plot_mz_decs <- function(mRList=NULL){

  
  quality <- data.frame(quality=mRList$metab_ann$mz %% 1)
 
  set_theme(theme_science())
  
  p <- ggplot(quality, aes(x=quality)) +
    geom_histogram(fill="#2C3E50", color="white", bins = 20) +
    ylab("frequency") +
    xlab("m/z decimal part")
 
  return(p)
}
