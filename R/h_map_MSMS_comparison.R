#' Heatmaps of errors between all-prais of peaks in feature and metabolite MS/MS spectra
#' @param mRList mRList object with metabolite_identification element
#' @param metab_id metabolite ID
#' @param feature_id feature ID
#' @param col color palette
#' @param ppm_error ppm_errors above this threshold are not shown
#' @param na_col color for NA values
#' @param RI_difference RI_difference above this threshold are not shown
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom pals brewer.purples
#' @importFrom grid gpar grid.text
#' @returns List of ComplexHeatmap objects
#' @export

h_map_MSMS_comparison <- function(mRList=NULL, metab_id=NULL, feature_id=NULL, col=NULL, na_col="black", ppm_error=20, RI_difference=20){
  
  if(is.null(col)){
    col <- brewer.purples(7)
  }
  
  cell_ppm <- function(j, i, x, y, w, h, fill) {
    if(is.na(X[i, j])){
      ppm_error_ij <- paste0(">", ppm_error)
      ppm_error_ij <- ""
      col_ij <- "black"
    }else{
      ppm_error_ij <- format(X[i, j], digits=1,  scientific = F)
      col_ij <- "red"
    }
    grid.text(ppm_error_ij, x, y, gp = gpar(fontsize = 9, col=col_ij, font=2))
  }
  
  id_assoc <- unique(mRList$metabolite_identification$associations[mRList$metabolite_identification$associations$ID == metab_id & mRList$metabolite_identification$associations$Feature_ID == feature_id, c("ID", "Feature_ID", "ID_peaks"), drop=F])
  
  if(nrow(id_assoc)==0){
    stop("Can't find associations between ", metab_id, " and ", feature_id, ".\n")
  }
  
  p_list <- vector("list", nrow(id_assoc))
  for(i in 1:nrow(id_assoc)){
    
    m_list <- mRList$metabolite_identification$MS_MS_info[[id_assoc$ID[i]]][[id_assoc$ID_peaks[i]]][[id_assoc$Feature_ID[i]]]
    
    if (is.null(m_list)) {
      stop(paste0("Considering the passed ", metab_id, " and ", feature_id, ", there is no suitable ms/ms information"))
    }
    
    X <- m_list$ppm_error
    
    X[X > ppm_error] <- NA
    
    hm <- Heatmap(X, cluster_rows = F, cluster_columns = F, row_labels = 1:nrow(X), column_labels = 1:ncol(X), name = "e_mz", col=col, cell_fun = cell_ppm, rect_gp = gpar(col = "white", lwd = 1, lty=2), row_names_side = "left", row_title = paste(id_assoc$ID[i], id_assoc$ID_peaks[i]), column_title = id_assoc$Feature_ID[i], column_title_side = "bottom", na_col = na_col)
    p_list[[i]] <- list(e_mz=hm)
    
    cell_RI <- function(j, i, x, y, w, h, fill) {
      if(is.na(X[i, j])){
        #RI_difference_ij <- "paste0(">", RI_difference)"
        RI_difference_ij <- ""
        col_ij <- "black"
      }else{
        RI_difference_ij <- format(X[i, j], digits=1,  scientific = F)
        col_ij <- "red"
      }
      grid.text(RI_difference_ij, x, y, gp = gpar(fontsize = 9, col=col_ij, font=2))
    }
    
    X <- m_list$RI_diff
    X[X > RI_difference] <- NA
    
    hm <- Heatmap(X, cluster_rows = F, cluster_columns = F, row_labels = 1:nrow(X), column_labels = 1:ncol(X), name = "e_I", col=col, cell_fun = cell_RI, rect_gp = gpar(col = "white", lwd = 1, lty=2), row_names_side = "left", row_title = paste(id_assoc$ID[i], id_assoc$ID_peaks[i]), column_title = id_assoc$Feature_ID[i], column_title_side = "bottom", na_col = na_col)
    p_list[[i]]$e_RI <- hm
    
  }
  
  return(p_list)
  
}
