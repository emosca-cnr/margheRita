#PCA function post processing
#Output score, loading plot (PC1 vs PC2) and scree plot
#@param df_post: dataframe of data was processed:imputation of missing value, pareto scaled




#PCA post processing directory
#dirout = paste(getwd(), "/PCA_Post", "/", sep = "")
#dir.create(dirout)

pca_postprocessing <- function(df_post,dirout) {
  dirout = paste(getwd(), "/PCA_Post", "/", sep = "")
  dir.create(dirout)
  pca.post <- prcomp(t(df_post[,4:ncol(df_post)]), scale = T, center = T)
  p.v.post = matrix(((pca.post$sdev ^ 2) / (sum(pca.post$sdev ^ 2))), ncol = 1) #varianza
  p.i.post = round(p.v.post * 100, 1) #percentuali di varianza spiegata dalle PC
  pwd.score.post = paste(getwd(), "/PCA_Post/PCA_Post_ScoreMatrix.csv", sep =
                           "")
  write.csv(pca.post$x, pwd.score.post)
  pwd.load.post = paste(getwd(), "/PCA_Post/PCA_Post_LoadingsMatrix.csv", sep =
                          "")
  write.csv(pca.post$rotation, pwd.load.post)
  pwd.pvar.post = paste(getwd(), "/PCA_Post/PCA_Post_Variance.csv", sep =
                          "")
  write.csv(p.i.post, pwd.pvar.post)
  Pvar.post = p.i.post
  #ora faccio i grafici di score, loading and scree plot plotting PC1 vs PC2
  scoreplot_post = paste(dirout, "Scoreplot_post.png", sep = "")
  png(
    scoreplot_post,
    width = 8,
    height = 8,
    units = "in",
    res = 300
  )
  plot(
    pca.post$x[, 1],
    pca.post$x[, 2],
    xlab = paste("PC1 (", p.i.post[, 1], "%)", sep = ""),
    ylab = paste("PC2 (", p.i.post[, 2], "%)", sep = ""),
    main = "Score plot post",
    col = metadata$Group,
    pch = 19
  )
  dev.off()
  loadingplot_post = paste(dirout, "Loadingplot_post.png", sep = "")
  png(
    loadingplot_post,
    width = 8,
    height = 8,
    units = "in",
    res = 300
  )
  plot(
    pca.post$rotation[, 1],
    pca.post$rotation[, 2],
    xlab = "Loading 1",
    ylab = "Loading 2",
    main = "Loading plot post",
    pch = 19
  )
  dev.off()
  scree_post = paste(dirout, "Screeplot_post.png", sep = "")
  png(
    scree_post,
    width = 8,
    height = 8,
    units = "in",
    res = 300
  )
  barplot(
    Pvar_post[, 1],
    xlab = "Principal Components",
    ylab = "Proportion of Variance explained",
    main = "Screeplot_post",
    ylim = c(0, 100)
  )
  dev.off()
}



