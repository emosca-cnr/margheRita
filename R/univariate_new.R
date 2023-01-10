#' Univariate analysis
#' @description User can choose between T-test, U-test, Anova test  or Kruskal-Wallis test.
#' @param mRList dataframe with data of samples that contains only biological replicates. Technical replicates have to be collapsed before using collapse_tech_rep function
#' @param test_method choice of the test to apply
#' @param paired TRUE if data are paired, FALSE for unpaired data
#' @param group_factor specify the groups to compare
#' @param dirout output directory
#' @importFrom utils write.csv
#' @importFrom stats na.omit t.test p.adjust wilcox.test anova aov kruskal.test
#' @export
#' @return mRList object with mRList$"testchosen" with univariate analysis
#' @examples
#' ##library(dataset.margheRita)
#' ##dataset(norm_pos)
#' mRList<-univariate(mRList, test_method="anova",paired="FALSE",group_factor="class")


univariate_new<- function(mRList, dirout="./", test_method=c("ttest","Utest", "anova","kruskal"), paired=c("FALSE", "TRUE"), contrast_samples="class"){

  paired<- match.arg(paired)
  test_method <- match.arg(test_method)

  dirout = paste(dirout, sep = "")
  dir.create(dirout)


  X_ann <- mRList$sample_ann
  X_ann <- na.omit(X_ann)
  X_ann$contrast_samples<-X_ann$class
  X_ann$contrast_samples[! X_ann$contrast_samples %in% contrast_samples] <- NA
  X_ann <- na.omit(X_ann)
  group_factor <- as.factor(X_ann$contrast_samples)
  X_data<-mRList$data[,rownames(X_ann)]

  if(test_method == "ttest"){

    if(length(levels(group_factor))>2){
      cat("groups:", levels(group_factor), "\n")
      stop("can not apply t test for more than 2 groups")
    }

    cat("t tests between", levels(group_factor), "\n")

    ans <- apply(X_data, 1, function(x) t.test(x ~ group_factor)$ "p.value"[1])
    ans<-as.data.frame(ans)
    colnames(ans)<-"pvalue"
    #ans <- lapply(ans, function(x) data.frame(estimate=x$estimate, stderr=x$stderr, t=x$statistic, p=x$p.value))
    #ans <- do.call(rbind, ans_a)
    ans$qvalue <- p.adjust(ans$pvalue, method ="BH")
    #colnames(ans)<-c("t-statistics","p-value","q-value")
    mRList$ttest <- ans
    #rownames(mRList$ttest)<-rownames(mRList$data)
    uni_t= paste(dirout, "/ttest.csv", sep ="")
    utils::write.csv(mRList$ttest, uni_t)
  }

  if(test_method == "Utest"){

    if(length(levels(group_factor))>2){
      cat("groups:", levels(group_factor), "\n")
      stop("can not apply Utest for more than 2 groups")
    }

    cat("U tests between", levels(group_factor), "\n")

    ans <- apply(X_data, 1, function(x)  wilcox.test(x ~ group_factor))
    ans <- lapply(ans, function(x) data.frame(t=x$statistic, p=x$p.value))
    ans <- do.call(rbind, ans)
    ans$q <- p.adjust(ans$p, method ="BH")
    colnames(ans)<-c("t-statistics","p-value","q-value")
    mRList$Utest <- ans
    rownames(mRList$Utest)<-rownames(mRList$data)
    uni_u= paste(dirout, "/Utest.csv", sep ="")
    utils::write.csv(mRList$Utest, uni_u)

  }



  if(test_method == "anova"){

    if(length(levels(group_factor))<3){
      cat("groups:", levels(group_factor), "\n")
      stop("can not apply anova test for less than 3 groups")
    }

    ans_t <- apply(X_data, 1, function(x)  anova(aov(x ~ group_factor))$"F value"[1])
    ans_p<- apply(X_data, 1, function(x)  anova(aov(x ~ group_factor)) $"Pr(>F)"[1])
    ans_t<-as.data.frame(ans_t)
    ans_p<-as.data.frame(ans_p)
    anova_res<-cbind(ans_t,ans_p)
    colnames(anova_res)<-c("t","p")
    anova_res$q <- p.adjust(anova_res$p, method ="fdr")
    anova_res<-cbind(anova_res$t,anova_res$p,anova_res$q)
    colnames(anova_res)<-c("t-statistic","p-value","q-value")
    anova_res<-as.data.frame(anova_res)
    mRList$anova <- anova_res
    rownames(mRList$anova)<-rownames(mRList$data)

    uni_anova= paste(dirout, "/anova.csv", sep ="")
    utils::write.csv(mRList$anova, uni_anova)
  }

  if(test_method == "kruskal"){

    if(length(levels(group_factor))<3){
      cat("groups:", levels(group_factor), "\n")
      stop("can not apply kruskal test for less than 3 groups")
    }

    #cat("U tests between", levels(group_factor), "\n")

    ans <- apply(X_data, 1, function(x)  kruskal.test(x ~ group_factor))
    ans <- lapply(ans, function(x) data.frame(t=x$statistic, p=x$p.value))
    ans <- do.call(rbind, ans)
    ans$q <- p.adjust(ans$p, method ="fdr")
    colnames(ans)<-c("t-statistics","p-value","q-value")
    mRList$kruskal <- ans
    rownames(mRList$kruskal)<-rownames(mRList$data)
    uni_kru= paste(dirout, "/kruskal_test.csv", sep ="")
    utils::write.csv(mRList$kruskal, uni_kru)
  }

  return(mRList)

}


