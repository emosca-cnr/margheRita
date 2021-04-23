#' @import limma
#' @export
#Remove batch effect
#@param remove batch effect from df (dataframe of data)

#library(limma)
#boxplot(df[4:NCOL(df)], main="Original")
#df <- removeBatchEffect(df[4:NCOL(df)], batch)
#boxplot(df[4:NCOL(df)], main="Batch corrected")
bath_effect <- function(df, columns){
  box_pre <- boxplot(df[columns[1]:columns[2]], main="Original")
  df <- removeBatchEffect(df[columns[1]:columns[2]], batch)
  box_post <- boxplot(df[4:ncol(df)], main="Batch corrected")
  box_plot <- box_pre +  box_post
  return(as.data.frame(df), box_plot)
}



