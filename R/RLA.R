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

RLA <- function(mRList, include_QC=FALSE, logged=FALSE, robust=TRUE, do_plot=FALSE, out_dir="./", ...){


  if(include_QC){
    mRList$data <- cbind(mRList$data, mRList$QC)
    mRList$sample_ann<- rbind(mRList$sample_ann, mRList$QC_ann)
  }


  ans <- mRList$data
  if(!logged){
    ans <- log2(mRList$data + 1)
  }

	if(robust){
		ans <- t(apply(ans, 1, function(y) y - stats::median(y)))
	}else{
		ans <- t(apply(ans, 1, function(y) y - mean(y)))
	}

	if(do_plot){
	  ## add col_by
	  dir.create(out_dir)
	  jpeg(paste0(out_dir, "/RLA.jpg"), width = 200, height = 100, res=300, units="mm")
	  par(mar=c(4, 4, 1, 1))
	  boxplot(ans, ..., ylab="x - <x>")
	  abline(h=0, lty=2)
	  dev.off()
	}

	return(ans)

}
