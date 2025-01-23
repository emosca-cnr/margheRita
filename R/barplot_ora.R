
#' Draw a bar plot from pathway analyses with ORA
#' @param pa_res pa_res obtained from pathway_analysis with type "ora".
#' @param p "pvalue" or "p.adjust".
#' @param filename Name of the .png file that will be created.
#' @export
#' @importFrom ggplot2 ggplot ggsave

barplot_ora <- function(pa_res, p = c("pvalue", "p.adjust"), filename = "pathway_analysis_ora") {
  
  p <- match.arg(p, c("pvalue", "p.adjust"))
  
  pa_res_ora_table <- pa_res[["res"]]@result
  
  the_bar_plot <- ggplot(data = pa_res_ora_table, aes(x = Description, y = Count, fill = !!sym(p))) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "red", high = "blue", limits = c(0, 1)) +
    coord_flip() + 
    ggtitle("Over Representation Analysis (ORA)") +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
          plot.title = element_text(hjust = 0.5))
  
  ggsave(filename = paste0(filename, ".png"), plot = the_bar_plot, device = "png")
  
}

