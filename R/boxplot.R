#'@boxplot
#'@color and names of samples depends on m_list$sample_ann$class
#' @param m_list margheRita m_list
#' @export
#' @importFrom graphics plot barplot
#' @importFrom grDevices dev.off png




boxplot<-function(m_list, dirout){

  dirout = paste(dirout, sep = "")
  dir.create(dirout)

  col_factor <- as.factor(m_list$sample_ann$class)
  col_pal <- rainbow(length(levels(col_factor)))

  pdf("Boxplot.pdf")
    for (i in 1:nrow(m_list$data)) {
      boxplot(
        m_list$data[i,] ~ as.factor(m_list$sample_ann$class) ,
        names =levels(as.factor(m_list$sample_ann$class)),
        #main = m_list$metabo_ann, #metabolites or Ms ID or annotation
        ylab="Relative Abundance", col= col_pal)
    }
dev.off()
}
