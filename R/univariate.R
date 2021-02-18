#Shapiro-Wilk's test is used for testing the normal distribution of metabolites lo faccio sul subsetting Samples
#after technical replicates collapse
#counting group
#univariate on df_sample_mean only on biological replicates
#@param df_samples_mean (dataframe of samples only biological replicates, technical replicates were collapsed before)

conteggio = nlevels(as.factor(metadata_samples$Group))

univariate_analysis <- function(df_samples_mean) {
  uni = c()
  test = c()
  for (i in 1:nrow(df_samples_mean, metadata_samples$Group)) {
    pvalueShapiro <- shapiro.test(df_samples_mean[i, ])
    if (pvalueShapiro$p.value < 0.05) {
      print("it is not normal distributed")
      if (conteggio == 2) {
        uni <-
          c(uni,
            wilcoxon.test(df_samples_mean[i, ] ~ metadata_samples$Group))
        test <- c(test, "wilcoxon")
      }
      else {
        uni <- c(uni,
                 kruskal.test(df_samples_mean[i, ] ~ metadata_samples$Group))
        test <- c(test, "kruskal")
      }
    }
    else{
      print("it is normal distributed")
      if (conteggio == 2) {
        uni <- c(uni,
                 t.test(df_samples_mean[i, ] ~ metadata_samples$Group))
        test <- c(test, "ttest")
      }
      else {
        uni <- c(uni, aov(df_samples_mean[i, ] ~ metadata_samples$Group))
        test <- c(test, "anova")
      }
    }
  }
  uni_corrected <- BH(uni)
  matrix_analysis <- cbind(df_samples_mean, uni, uni_corrected, test)
  return(matrix_analysis)
  return(uni)
  return(uni_corrected)
  write.csv(matrix_analysis, "univariate_analysis_samples.csv")
}

univariate_analysis(df_samples_mean)


#heatmap on samples (no QC) using mean matrix in whch column are samples and rows metabolites

library("RColorBrewer")
library(pheatmap)

pdf("heatmap.pdf")

col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
annotation_column <- metadata_samples[, 1:ncol(metadata_samples)]
mycolors_s <-
  brewer.pal(n = conteggio, name = "Dark2")
names(mycolors_s) = levels(metadata_samples$Group)
ann_colors = list(samples = mycolors_s)
crp <- colorRampPalette(c('blue', 'white', 'red'))
colors = nrow(df_samples_mean)
pheatmap(
  df_samples_mean,
  annotation_col = annotation_column,
  annotation_colors = ann_colors,
  col = colors,
  cellwidth = 20,
  cellheight = 1.60,
  show_rownames = FALSE,
  labels_row = row.names(df_samples_mean),
  legend = TRUE
)

dev.off()



#Boxplot only for p-values<0.05 Boxplot for each metabolite (row) or uni_corrected?

pdf("Boxplot.pdf")
for (i in 1:nrow(df_samples_mean)) {
  if (uni[i] < 0.05) {
    boxplot(
      df_samples_mean[i, ] ~ metadata_samples$Group,
      col = metadata_samples$Group,
      names = row.names(samples_names),
      main = row.names(df_samples_mean)
    )
  }
}
dev.off()
