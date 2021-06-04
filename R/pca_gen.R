#' PCA general
#' Scree plot, score and loading plot (PC1 vs PC2)
#' @param m_list margheRita m_list
#' @export
#' @importFrom graphics plot barplot
#' @importFrom grDevices dev.off png
#' @importFrom utils write.csv


#PCA function general to use in different points

pca_gen <- function(m_list, dirout, col_by="class", scale=TRUE, include_QC=TRUE) {
  dirout = paste(dirout, sep = "")
  dir.create(dirout)
  if (scale) {
    #m_list$data<-apply(m_list$data,1,function(x) pareto(x))
    m_list <- pareto(m_list) #ettore: pareto requires m_list
  }
  if (include_QC){
    data_ <- cbind(m_list$data, m_list$QC)
  }
  pca <- prcomp(t(data_), scale = F, center = F)
  p.v.= matrix(((pca$sdev ^ 2) / (sum(pca$sdev ^ 2))), ncol = 1) #varianza
  p.i. = round(p.v.* 100, 1) #percentuali di varianza spiegata dalle PC
  pwd.score= paste(dirout, "/ScoreMatrix.csv", sep ="")
  utils::write.csv(pca$x, pwd.score)
  pwd.load. = paste(dirout ,"/LoadingsMatrix.csv", sep= "")
  utils::write.csv(pca$rotation, pwd.load.)
  pwd.pvar.= paste(dirout, "/Variance.csv", sep = "")
  utils::write.csv(p.i., pwd.pvar.)
  Pvar. = p.i.
  #ora faccio i grafici di score, loading and scree plot plotting PC1 vs PC2
  scoreplot = paste(dirout, "/Scoreplot.png", sep = "")
  grDevices::png(
    scoreplot,
    width = 8,
    height = 8,
    units = "in",
    res = 300
  )
  col_factor <- as.factor(m_list$sample_ann[, col_by])
  col_pal <- rainbow(length(levels(col_factor)))
  graphics::plot(
    pca$x[, 1],
    pca$x[, 2],
    xlab = paste("PC1 (", p.i.[1], "%)", sep = ""),
    ylab = paste("PC2 (", p.i.[2], "%)", sep = ""),
    main = "Score plot",
        col = col_pal[as.numeric(col_factor)],
        pch = 19
  )
  legend("bottomright", legend = levels(col_factor), col = col_pal, pch=16, cex=0.5)

  grDevices::dev.off()
  loadingplot = paste(dirout, "/Loadingplot.png", sep = "")
  grDevices::png(
    loadingplot,
    width = 8,
    height = 8,
    units = "in",
    res = 300
  )
  graphics::plot(
    pca$rotation[, 1],
    pca$rotation[, 2],
    xlab = "Loading 1",
    ylab = "Loading 2",
    main = "Loading plot",
    pch = 19
  )
  grDevices::dev.off()
  scree = paste(dirout, "/Screeplot.png", sep = "")
  grDevices::png(
    scree,
    width = 8,
    height = 8,
    units = "in",
    res = 300
  )
  graphics::barplot(
    Pvar.[, 1],
    xlab = "Principal Components",
    ylab = "Proportion of Variance explained",
    main = "Screeplot",
    ylim = c(0, 100)
  )
  grDevices::dev.off()
  return(m_list)
}
