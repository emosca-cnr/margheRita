#' user has to Choose between T-test, U-test, Anova or Kruskal
#' @param m_list$data (dataframe of samples only biological replicates, technical replicates were collapsed before)
#' @importFrom utils write.csv
#' @export


univariate <- function(m_list, dirout="./", test_method=c("ttest","anova","Utest","kruskal"),paired=c("FALSE","TRUE"), group_factor="class"){

  paired<- match.arg(paired)
  test_method <- match.arg(test_method)

  dirout = paste(dirout, sep = "")
  dir.create(dirout)

  X_ann <- m_list$sample_ann
  group_factor <- as.factor(X_ann[, group_factor])

  data<-t(m_list$data)


  uni = c()
  test = c()
  for (i in 1:nrow(m_list$data)){
    if (test_method=="ttest"){
      uni <- c(uni,t.test(m_list$data[i, ] ~ group_factor, data=m_list$data, paired=paired))
      test <- c(test, "ttest")
    }
    if (test_method=="anova"){
      uni <- c(uni, anova(m_list$data[i, ] ~ group_factor))
      test <- c(test, "anova")
    }
    if (test_method=="Utest"){
      uni <-c(uni, data=m_list$data,wilcox.test(m_list$data[i, ] ~ group_factor, paired=paired))
      test <- c(test, "wilcoxon")
    }
    if(test_method=="kruskal"){
      uni <- c(uni,kruskal.test(m_list$data[i, ] ~ group_factor,data=m_list$data, paired=paired))
      test <- c(test, "kruskal")
    }
    uni_corrected=c()
    uni_corrected <- p.adjust(uni, method ="BH")
    m_list$uni <- cbind(m_list$data$ID, uni, uni_corrected)
    #m_list$data<-cbind(mlist$data, uni_corrected)
    return(m_list)
    return(uni_corrected)
    #write.csv(uni_corrected, "univariate_corrected.csv")
  }
}
