#' Relative Log Abudance
#' @param mRList mRList object
#' @param include_QC whether to include or not the QC samples
#' @param logged are the input data on log-scale or not?
#' @param robust whether to use the median or not
#' @param do_plot whether to plot or not
#' @param out_dir output directory
#' @param ... further arguments to boxplot function
#' @export
#' @author Ettore Mosca (CNR-ITB)
#' @importFrom stats median
#' @importFrom grDevices jpeg rainbow
#' @importFrom graphics par abline
#' @importFrom pals alphabet2

RLA <- function(mRList=NULL, include_QC=FALSE, logged=FALSE, robust=TRUE, do_plot=FALSE, out_file="RLA.jpg", pal=NULL, col_by="class", ...){
  
  
  if(include_QC){
    X_data <- cbind(mRList$data, mRList$QC)
    X_sample_ann <- rbind(mRList$sample_ann, mRList$QC_ann)
  }else{
    X_data <- mRList$data
    X_sample_ann <- mRList$sample_ann
  }
  
  col_factor <- as.factor(X_sample_ann[, col_by])
  n <- length(levels(col_factor))
  
  if(is.null(pal) | length(pal) != n){
    pal <- pals::alphabet(n)
  }
  colors <- pal[as.numeric(col_factor)]
  
  if(!logged){
    X_data <- log2(X_data + 1)
  }
  
  if(robust){
    X_data <- t(apply(X_data, 1, function(y) y - stats::median(y)))
  }else{
    X_data <- t(apply(X_data, 1, function(y) y - mean(y)))
  }
  
  if(do_plot){
    ## add col_by
    jpeg(out_file, width = 200, height = 100, res=300, units="mm")
    par(mar=c(4, 4, 1, 1))
    boxplot(X_data, ..., ylab="x - <x>", main="Relative log Abudance", col = colors)
    abline(h=0, lty=2)
    
    legend("bottomright", legend = levels(col_factor), pch=15, col = pal, cex=0.5, horiz = T)
    dev.off()
  }
  
  mRList$RLA <- X_data
  
  return(mRList)
  
}
