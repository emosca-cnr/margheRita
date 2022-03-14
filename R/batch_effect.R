#' Remove batch effect
#' @importFrom limma removeBatchEffect
#' @export


batch_effect <- function(mRList){

  mRList$data <- limma::removeBatchEffect(mRList$data, mRList$sample_ann$batch)

  return(mRList)

}



