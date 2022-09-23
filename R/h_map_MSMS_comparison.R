#' Vh_map_MSMS_comparison
#' @param matrices_list list of matrices with ppm error, RI difference and flags
#' @param out_prefix prefix for output files
#' @importFrom ComplexHeatmap
#' @import viridis

h_map_MSMS_comparison <- function(matrices_list=NULL, out_prefix="MSMS_comparison", col=NULL, ppm_error=10, RI_difference=10){
  
  
  if(is.null(col)){
    col <- rev(viridis::viridis(7))
  }
  
  cell_ppm <- function(j, i, x, y, w, h, fill) {
    if(X[i, j] > ppm_error){
      ppm_error_ij <- paste0(">", ppm_error)
      ppm_error_ij <- ""
      col_ij <- "black"
    }else{
      ppm_error_ij <- format(X[i, j], digits=0,  scientific = F)
      col_ij <- "red"
    }
    grid.text(ppm_error_ij, x, y, gp = gpar(fontsize = 9, col=col_ij, font=2))
  }
  
  X <- matrices_list$ppm_error
  X[X > ppm_error] <- ppm_error+1
  
  jpeg(filename = paste0(out_prefix, "_ppm_error.jpg"), width = 200, height = 200, res=300, units = "mm")
  hm <- Heatmap(X, cluster_rows = F, cluster_columns = F, row_labels = 1:nrow(X), column_labels = 1:ncol(X), name = paste0("ppm\nerror\n<", ppm_error), col=col, cell_fun = cell_ppm, rect_gp = gpar(col = "white", lwd = 1, lty=2), row_names_side = "left", row_title = "library", column_title = "sample", column_title_side = "bottom")
  draw(hm)
  dev.off()
  
  cell_RI <- function(j, i, x, y, w, h, fill) {
    if(X[i, j] > RI_difference){
      #RI_difference_ij <- "paste0(">", RI_difference)"
      RI_difference_ij <- ""
      col_ij <- "black"
    }else{
      RI_difference_ij <- format(X[i, j], digits=1,  scientific = F)
      col_ij <- "red"
    }
    grid.text(RI_difference_ij, x, y, gp = gpar(fontsize = 9, col=col_ij, font=2))
  }
  
  X <- matrices_list$RI_diff
  X[X > RI_difference] <- RI_difference+1
  
  jpeg(filename = paste0(out_prefix, "_RI_diff.jpg"), width = 200, height = 200, res=300, units = "mm")
  hm <- Heatmap(X, cluster_rows = F, cluster_columns = F, row_labels = 1:nrow(X), column_labels = 1:ncol(X), name = paste0("RI\ndiff\n<", RI_difference), col=col, cell_fun = cell_RI, rect_gp = gpar(col = "white", lwd = 1, lty=2), row_names_side = "left", row_title = "library", column_title = "sample", column_title_side = "bottom")
  draw(hm)
  dev.off()
  
  
}