#' Relative Log Abudance
#' @param mRList mRList object
#' @param include_QC whether to include or not the QC samples
#' @param logged are the input data on log-scale or not?
#' @param robust whether to use the median or not
#' @param pal color palette
#' @param col_by column of sample_ann to color
#' @param ... further arguments to boxplot function
#' @export
#' @importFrom stats median
#' @importFrom pals alphabet2
#' @importFrom graphics abline
#' @return mRList object with element "RLA"

RLA <- function(mRList=NULL, include_QC=FALSE, logged=FALSE, robust=TRUE, pal=NULL, col_by="class", ...){


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
    pal <- alphabet2(n)
  }
  colors <- pal[as.numeric(col_factor)]

  if(!logged){
    X_data <- log2(X_data + 1)
  }

  if(robust){
    X_data <- t(apply(X_data, 1, function(y) y - median(y)))
  }else{
    X_data <- t(apply(X_data, 1, function(y) y - mean(y)))
  }

  ## add col_by
  boxplot(X_data, ..., ylab="x - <x>", main="Relative log abundance", col = colors)
  abline(h=0, lty=2)
  legend("bottomright", legend = levels(col_factor), pch=15, col = pal, cex=0.5, horiz = T)


  mRList$RLA <- X_data

  return(mRList)

}
