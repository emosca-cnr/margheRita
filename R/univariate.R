#' Univariate analysis
#' @description It applies T test, U test, Anova or Kruskal-Wallis test over dataset features.
#' @param mRList mRList object
#' @param test_method any between "ttest", "Utest", "anova" or "kruskal".
#' @param exp.factor column of mRList$sample_ann that defines the level for each sample
#' @param exp.levels the levels to be considered in the column specified by exp.factor
#' @param dirout output directory
#' @importFrom utils write.csv
#' @importFrom stats t.test p.adjust wilcox.test anova aov kruskal.test TukeyHSD
#' @export
#' @return mRList object with mRList$"testchosen" with univariate analysis

univariate <- function(mRList=NULL, dirout="./", test_method=c("ttest","Utest", "anova","kruskal"), exp.levels=NULL, exp.factor="class"){
  
  test_method <- match.arg(test_method)
  
  dir.create(dirout)
  
  cat("exp.levels: ", exp.levels, "\n")
  cat("exp.factor: ", exp.factor, "\n")
  
  idx_samples <- mRList$sample_ann[, exp.factor] %in% exp.levels
  exp_design <- data.frame(id=rownames(mRList$sample_ann)[idx_samples], level=factor(mRList$sample_ann[idx_samples, exp.factor]))
  cat("selecxted samples:\n")
  print(exp_design)
  
  X_data <- mRList$data[, match(exp_design$id, colnames(mRList$data))]
  
  if(test_method == "ttest"){
    
    if(length(levels(exp_design$level))>2){
      cat("groups:", levels(exp_design$level), "\n")
      stop("can not apply t test for more than 2 groups")
    }
    
    cat("t tests between", levels(exp_design$level), "\n")
    
    ans <- apply(X_data, 1, function(x) t.test(x ~ exp_design$level))
    ans <- lapply(ans, function(x) data.frame(t=x$statistic, p=x$p.value))
    ans <- do.call(rbind, ans)
    ans$q <- p.adjust(ans$p, method ="fdr")
    mRList$ttest <- ans
    write.csv(mRList$ttest, file=file.path(dirout, "t_test.csv"))
  }
  
  if(test_method == "Utest"){
    
    if(length(levels(exp_design$level))>2){
      cat("groups:", levels(exp_design$level), "\n")
      stop("can not apply Utest for more than 2 groups")
    }
    
    cat("U tests between", levels(exp_design$level), "\n")
    
    ans <- apply(X_data, 1, function(x)  wilcox.test(x ~ exp_design$level))
    ans <- lapply(ans, function(x) data.frame(U=x$statistic, p=x$p.value))
    ans <- do.call(rbind, ans)
    ans$q <- p.adjust(ans$p, method ="fdr")
    mRList$Utest <- ans
    write.csv(mRList$Utest, file = file.path(dirout, "U_test.csv"))
    
  }
  
  
  if(test_method == "anova"){
    
    if(length(levels(exp_design$level))<3){
      cat("groups:", levels(exp_design$level), "\n")
      stop("can not apply anova test for less than 3 groups")
    }
    
    cat("Performing ANOVA between", levels(exp_design$level),"\n")
    aov_res <- apply(X_data, 1, function(x)  aov(x ~ exp_design$level))
    anova_res <- lapply(aov_res, anova)
    anova_res <- do.call(rbind, lapply(anova_res, function(x) data.frame(F=x$`F value`[1], p=x$`Pr(>F)`[1])))
    tukey_res <- lapply(aov_res, TukeyHSD)
    
    anova_res$q <- p.adjust(anova_res$p, method ="fdr")
    
    tukey_p_adj <- do.call(rbind, lapply(tukey_res, function(x) x[[1]][, "p adj"]))
    anova_res <- merge(anova_res, tukey_p_adj, by=0, all=T, sort=F)
    rownames(anova_res) <- anova_res[, 1]
    anova_res[, 1] <- NULL
    
    mRList$anova <- anova_res
    mRList$anova <- mRList$anova[match(mRList$metab_ann$Feature_ID, rownames(mRList$anova)), ]
    
    write.csv(mRList$anova, file=file.path(dirout, "anova_test.csv"))
  }
  
  if(test_method == "kruskal"){
    
    if(length(levels(exp_design$level))<3){
      cat("groups:", levels(exp_design$level), "\n")
      stop("can not apply kruskal test for less than 3 groups")
    }
    
    ans <- apply(X_data, 1, function(x)  kruskal.test(x ~ exp_design$level))
    ans <- lapply(ans, function(x) data.frame(H=x$statistic, p=x$p.value))
    ans <- do.call(rbind, ans)
    ans$q <- p.adjust(ans$p, method ="fdr")
    mRList$kruskal <- as.data.frame( ans)
    write.csv(mRList$kruskal, file = file.path(dirout, "kruskal_test.csv"))
  }
  
  return(mRList)
  
}


