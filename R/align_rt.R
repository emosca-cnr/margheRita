#'Align metabolites from different batches on the basis of retention times
#
#' @param X data.frame with the following columns: PeakName, RT
#' @param rt_err maximum time difference between retention times of the same metabolite
#' @param rt_class internal parameter set by default to rt_err/2

align_rt <- function(X=NULL, rt_err=1, rt_class=NULL){
  
  
  map_retention_times <- function(X, rt_class=NULL){
    
    #1: peak name
    #2: RT
    
    X <- X[order(X$Peak.Name, X$RT), ]
    X$RT_class <- 0
    class <- 0
    curr_ret_time <- X$RT[1]
    for(i in 2:(nrow(X))){
      
      if(X$Peak.Name[i] != X$Peak.Name[i-1]){ #Peak Names
        class <- 0
      }else{
        if(abs(X$RT[i] - X$RT[i-1]) > rt_class){ #RT
          class <- class+1
        }
      }
      X$RT_class[i] <- class
    }
    
    class_mean <- lapply(split(X, X$Peak.Name), function(imetab) tapply(imetab$RT, imetab$RT_class, mean))
    class_mean <- data.frame(Peak.Name=rep(names(class_mean), times=unlist(lapply(class_mean, length))), RT_class=unlist(lapply(class_mean, names)), RT_mean=unlist(class_mean), stringsAsFactors = F)
    
    ans <- merge(X, class_mean, by=c("Peak.Name", "RT_class"), sort=F)
    
    return(ans)
    
  }
  
  

  if(is.null(rt_class)){
    rt_class <- rt_err #/ 2
  }
  
  cat("rt_err", rt_err, "\n")
  cat("rt_class", rt_class, "\n")
  
  #unique peaks and retention times
  #pair retantion time
  ans <- data.frame(X[, 1:2], n_mes=rowSums(sign(X[, -c(1:2)]), na.rm = T), stringsAsFactors = F)
  
  #force colnames
  colnames(ans)[1:2] <- c("Peak.Name", "RT")
  
  #mapr retention times within eps minutes
  ans <- merge(ans, map_retention_times(ans[, c(1:2)], rt_class = rt_class), by=c("Peak.Name", "RT"))
  
  #check classes for retention times greater than rt_err 
  temp <- unlist(lapply(split(ans, ans$Peak.Name), function(x) tapply(x[, 2], x$RT_mean, function(y) ifelse((max(y) - min(y)) > rt_err, 1, 0))))
  if(any(temp>0)){
    cat("metabolites assigned to the same class beyond", rt_err, ". Try rt_class <", rt_class, "\n")
    print(temp[temp>0])
  }
  
  #total number of measurements in each rt class
  temp <- lapply(split(ans, ans[, 1]), function(x) tapply(x$n_mes, x$RT_class, sum))
  temp_perc <- lapply(temp, function(x) x/sum(x)) #fraction of measurements of each class
  temp <- data.frame(Peak.Name=rep(names(temp), times=unlist(lapply(temp, length))), RT_class=unlist(lapply(temp, names)), tot_n_mes=unlist(temp), perc_tot=unlist(temp_perc), stringsAsFactors = F)
  ans <- merge(ans, temp, by=c("Peak.Name", "RT_class"), sort=F)
  ans$perc_tot[is.nan(ans$perc_tot)] <- 0
  ans
  
}
