#' PCA after normalizazion
#' Scree plot, score and loading plot (PC1 vs PC2)
#' @param df_norm: dataframe after normalization
#' @export
#' @importFrom graphics plot barplot
#' @importFrom grDevices dev.off png
#' @importFrom utils write.csv

#PCA after normalization directory
#dirout = paste(getwd(), "/PCA_Norm/", sep = "")
#dir.create(dirout)

#PCA function after norm and plots

pca_norm <- function(df_norm,dirout) {
  dirout = paste(getwd(), "/PCA_Norm/", sep = "")
  dir.create(dirout)
  pca.norm <- prcomp(t(df_norm[,4:ncol(df_norm)]), scale = T, center = T)
  p.v.norm = matrix(((pca.norm$sdev ^ 2) / (sum(pca.norm$sdev ^ 2))), ncol = 1) #varianza
  p.i.norm = round(p.v.norm * 100, 1) #percentuali di varianza spiegata dalle PC
  pwd.score.norm = paste(getwd(), "/PCA_Norm/PCA_Norm_ScoreMatrix.csv", sep ="")
  utils::write.csv(pca.norm$x, pwd.score.norm)
  pwd.load.norm = paste(getwd(), "/PCA_Norm/PCA_Norm_LoadingsMatrix.csv", sep= "")
  utils::write.csv(pca.norm$rotation, pwd.load.norm)
  pwd.pvar.norm = paste(getwd(), "/PCA_Norm/PCA_Pre_Variance.csv", sep =
                          "")
  utils::write.csv(p.i.norm, pwd.pvar.norm)
  Pvar.norm = p.i.norm
  #ora faccio i grafici di score, loading and scree plot plotting PC1 vs PC2
  scoreplot_norm = paste(dirout, "scoreplot_norm.png", sep = "")
  grDevices::png(
    scoreplot_norm,
    width = 8,
    height = 8,
    units = "in",
    res = 300
  )
  graphics::plot(
    pca.norm$x[, 1],
    pca.norm$x[, 2],
    xlab = paste("PC1 (", p.i.norm[1], "%)", sep = ""),
    ylab = paste("PC2 (", p.i.norm[2], "%)", sep = ""),
    main = "Score plot norm",
    col = metadata$Group,
    pch = 19
  )
  grDevices::dev.off()
  loadingplot_norm = paste(dirout, "loadingplot_norm.png", sep = "")
  grDevices::png(
    loadingplot_norm,
    width = 8,
    height = 8,
    units = "in",
    res = 300
  )
  graphics::plot(
    pca.norm$rotation[, 1],
    pca.norm$rotation[, 2],
    xlab = "Loading 1",
    ylab = "Loading 2",
    main = "Loading plot norm",
    pch = 19
  )
  grDevices::dev.off()
  scree_norm = paste(dirout, "Screeplot_norm.png", sep = "")
  grDevices::png(
    scree_norm,
    width = 8,
    height = 8,
    units = "in",
    res = 300
  )
  graphics::barplot(
    Pvar.norm[, 1],
    xlab = "Principal Components",
    ylab = "Proportion of Variance explained",
    main = "Screeplot_norm",
    ylim = c(0, 100)
  )
  grDevices::dev.off()
  return(pca.norm)
}


