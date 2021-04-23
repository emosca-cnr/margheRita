#' @import limma
#' @export
#Remove batch effect
#@param remove batch effect from df (dataframe of data)

#library(limma)
#boxplot(df[4:NCOL(df)], main="Original")
#df <- removeBatchEffect(df[4:NCOL(df)], batch)
#boxplot(df[4:NCOL(df)], main="Batch corrected")
batch_effect <- function(m_list){

  #box_pre <- boxplot(list$df[, -c(1:3)], main="Original")
  m_list$df <- limma::removeBatchEffect(m_list$df[, -c(1:3)], m_list$metadata$batch)

  #box_post <- boxplot(df[, -c(1:3)], main="Batch corrected")
  #box_plot <- box_pre +  box_post

  return(m_list)

}



