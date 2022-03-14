#' user has to Choose between T-test, U-test, Anova or Kruskal
#' @param mRList (dataframe of samples only biological replicates, technical replicates were collapsed before)
#' @importFrom utils write.csv
#' @export


univariate <- function(mRList, dirout="./", test_method=c("ttest", "anova", "Utest", "kruskal"), paired=c("FALSE", "TRUE"), group_factor="class"){

  paired<- match.arg(paired)
  test_method <- match.arg(test_method)

  dirout = paste(dirout, sep = "")
  dir.create(dirout)

  X_ann <- mRList$sample_ann
  group_factor <- as.factor(X_ann[, group_factor])

  data<-t(mRList$data)


  uni = c()
  test = c()
  for (i in 1:nrow(mRList$data)){
    if (test_method=="ttest"){
      uni <- c(uni,t.test(mRList$data[i, ] ~ group_factor, data=mRList$data, paired=paired))
      test <- c(test, "ttest")
    }
    if (test_method=="anova"){
      uni <- c(uni, anova(mRList$data[i, ] ~ group_factor))
      test <- c(test, "anova")
    }
    if (test_method=="Utest"){
      uni <-c(uni, data=mRList$data,wilcox.test(mRList$data[i, ] ~ group_factor, paired=paired))
      test <- c(test, "wilcoxon")
    }
    if(test_method=="kruskal"){
      uni <- c(uni,kruskal.test(mRList$data[i, ] ~ group_factor,data=mRList$data, paired=paired))
      test <- c(test, "kruskal")
    }
    uni_corrected <- p.adjust(uni, method ="BH")
    mRList$uni <- cbind(mRList$data$ID, uni, uni_corrected)
    #mRList$data<-cbind(mlist$data, uni_corrected)
     #write.csv(uni_corrected, "univariate_corrected.csv")
  }

  return(mRList)

}
