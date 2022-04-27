#' user has to Choose between T-test, U-test, Anova or Kruskal
#' @param mRList (dataframe of samples only biological replicates, technical replicates were collapsed before)
#' @importFrom utils write.csv
#' @export


univariateEM <- function(mRList, dirout="./", test_method=c("ttest", "anova", "Utest", "kruskal"), paired=c("FALSE", "TRUE"), group_factor="class"){

  paired<- match.arg(paired)
  test_method <- match.arg(test_method)

  dirout = paste(dirout, sep = "")
  dir.create(dirout)

  X_ann <- mRList$sample_ann
  X_ann <- na.omit(X_ann)
  group_factor <- as.factor(X_ann[, group_factor])

  X_data <- mRList$data[, colnames(mRList$data) %in% rownames(X_ann)]

  if(test_method == "ttest"){
    
    if(length(levels(group_factor))>2){
      cat("groups:", levels(group_factor), "\n")
      stop("can not apply t test for more than 2 groups")
    }

    cat("t tests between", levels(group_factor), "\n")
    
    ans <- apply(X_data, 1, function(x) t.test(x ~ group_factor))
    ans <- lapply(ans, function(x) data.frame(t=x$statistic, p=x$p.value))
    ans <- do.call(rbind, ans)
    ans$q <- p.adjust(ans$p, method ="BH")
    
    mRList$ttest <- ans
  }

  
  return(mRList)

}
