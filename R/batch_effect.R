#' Removes batch effects
#' 
#' This function applies the procedure implemented in "removeBatchEffect" from limma package.
#' @param mRList mRList object
#' @importFrom limma removeBatchEffect
#' @export
#' @return mRList object with updated mRList$data.


batch_effect <- function(mRList){

  mRList$data <- limma::removeBatchEffect(mRList$data, mRList$sample_ann$batch)

  return(mRList)

}



