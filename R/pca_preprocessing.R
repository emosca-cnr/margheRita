#PCA function pre processing
#Output score, loading plot (PC1 vs PC2) and scree plot
#@param df: dataframe of data in which the first three column are MS DIAL ID, RT and m/z, sample in colum and features (metabolites) in rows


#PCA preprocessing directory
#dirout = paste(getwd(), "/PCA_Pre/", sep = "")
#dir.create(dirout)

#PCA function before processing and plots

pca_preprocessing <- function(df,dirout) {
  dirout = paste(getwd(), "/PCA_Pre/", sep = "")
  dir.create(dirout)
  pca.pre <- prcomp(t(df[, 4:ncol(df)]), scale = T, center = T)
  p.v.pre = matrix(((pca.pre$sdev ^ 2) / (sum(pca.pre$sdev ^ 2))), ncol = 1) #varianza
  p.i.pre = round(p.v.pre * 100, 1) #percentuali di varianza spiegata dalle PC
  pwd.score.pre = paste(getwd(), "/PCA_Pre/PCA_Pre_ScoreMatrix.csv", sep =
                          "")
  write.csv(pca.pre$x, pwd.score.pre)
  pwd.load.pre = paste(getwd(), "/PCA_Pre/PCA_Pre_LoadingsMatrix.csv", sep =
                         "")
  write.csv(pca.pre$rotation, pwd.load.pre)
  pwd.pvar.pre = paste(getwd(), "/PCA_Pre/PCA_Pre_Variance.csv", sep = "")
  write.csv(p.i.pre, pwd.pvar.pre)
  Pvar.pre = p.i.pre
  #ora faccio i grafici di score, loading and scree plot plotting PC1 vs PC2
  scoreplot_pre = paste(dirout, "Scoreplot_pre.png", sep = "")
  png(
    scoreplot_pre,
    width = 8,
    height = 8,
    units = "in",
    res = 300
  )
  plot(
    pca.pre$x[, 1],
    pca.pre$x[, 2],
    xlab = paste("PC1 (", p.i.pre[1], "%)", sep = ""),
    ylab = paste("PC2 (", p.i.pre[2], "%)", sep = ""),
    main = "Score plot pre",
    col = metadata$Group,
    pch = 19
  )
  dev.off()
  loadingplot_pre = paste(dirout, "Loadingplot_pre.png", sep = "")
  png(
    loadingplot_pre,
    width = 8,
    height = 8,
    units = "in",
    res = 300
  )
  plot(
    pca.pre$rotation[, 1],
    pca.pre$rotation[, 2],
    xlab = "Loading 1",
    ylab = "Loading 2",
    main = "Loading plot pre",
    pch = 19
  )
  dev.off()
  scree_pre = paste(dirout, "Screeplot_pre.png", sep = "")
  png(
    scree_pre,
    width = 8,
    height = 8,
    units = "in",
    res = 300
  )
  barplot(
    Pvar.pre[, 1],
    xlab = "Principal Components",
    ylab = "Proportion of Variance explained",
    main = "Screeplot_pre",
    ylim = c(0, 100)
  )
  dev.off()
}

