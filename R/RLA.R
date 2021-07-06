#' Relative Log Abudance
#' @param m_list margheRita list
#' @export
#' @author Ettore Mosca (CNR-ITB)
#' @importFrom stats median

RLA <- function(m_list, include_QC=FALSE, logged=FALSE, robust=TRUE, do_plot=FALSE, out_dir="./", ...){


  if(include_QC){
    m_list$data <- cbind(m_list$data, m_list$QC)
    m_list$sample_ann<- rbind(m_list$sample_ann, m_list$QC_ann)
  }


  ans <- m_list$data
  if(!logged){
    ans <- log2(m_list$data + 1)
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
	  boxplot(ans, ..., ylab="x - <x>")
	  dev.off()
	}

	return(ans)

}
