#' Remove batch effect
#' @import limma
#' @import mRList
#' @export


batch_effect <- function(mRList){

  mRList$data <- limma::removeBatchEffect(mRList$data, mRList$sample_ann$batch)

  return(mRList)

}



