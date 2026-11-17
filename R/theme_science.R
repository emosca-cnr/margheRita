#' A ggplot2 theme for good looking plots in science
#' @importFrom ggplot2 theme_bw theme element_blank element_text element_rect
#' @returns NULL
#' @export
theme_science <- function(){ 
	theme_bw() + #based on bw
		theme(
			panel.grid = element_blank(), #remove background grid
			text=element_text(color="black"), #black text
			axis.text=element_text(color="black"), #black text in axis tick labels
			panel.border = element_rect(fill=NA, colour = "black", linewidth=1), #increase plot box border size
			legend.text = element_text(size = 11),   # Legend labels
			legend.title = element_text(size = 11, face = "bold")
		)
}


