#' cluster_retention_times
#' @importFrom utils tail
#' @export

cluster_retention_times <- function(X, rt_class=NULL){

  #1: peak name
  #2: RT

  ans <- X
  colnames(ans) <- c("Peak.Name", "RT")
  ans$RT_mean <- ans$RT
  #order peaks by id and ascending RT
  ans <- unique(ans[order(ans$Peak.Name, ans$RT_mean), ])

  peaks <- unique(ans$Peak.Name)
  n_classes <- nrow(ans)
  n_classes_curr <- 0
  while(n_classes_curr != n_classes){ #STOP criteria: when no further aggregation take place

    #update the number of different classes
    n_classes <- nrow(unique(ans[, c("Peak.Name", "RT_mean")]))
    cat("n_classes, n_classes_curr")
    cat(n_classes, n_classes_curr, "\n")

    #RT_class to 0
    ans$RT_class <- 0
    #define peak name RT_mean pair
    ans$Peak.Name_RT_mean <- apply(ans[, c("Peak.Name", "RT_mean")], 1, function(irow) paste0(irow, collapse = "_"))
    ans$Peak.Name_RT_mean <- gsub(" ", "", ans$Peak.Name_RT_mean)

    #cycle through peaks
    for(i in 1:length(peaks)){
      #cat(i, "\n")

      #indices of the peak RT with different RT mean
      curr_idx <- which(ans$Peak.Name==peaks[i] & !duplicated(ans$Peak.Name_RT_mean))

      delta <- c(1000, 1000)

      if(length(curr_idx) > 1){

        #from the first to the penultimate
        for(j in 1:(length(curr_idx)-1)){

          #RT difference
          delta[2] <- ans$RT_mean[curr_idx[j+1]] - ans$RT_mean[curr_idx[j]]

          #if all delta are greater than limit, current assign the 0 to current metabolite
          if(all(delta > rt_class)){
            ans$RT_class[curr_idx[j]] <- 0
          }

          #if any delta is within permitted time
          if(any(delta < rt_class)){
            # if the forward is lower than backword, assign N, else P
            if(delta[2] < delta[1]){
              ans$RT_class[curr_idx[j]] <- "N"
            }else{
              ans$RT_class[curr_idx[j]] <- "P"
            }
          }
          delta[1] <- delta[2] #update delta for next window
        }

        #last element
        if(delta[2] > rt_class){
          ans$RT_class[curr_idx[length(curr_idx)]] <- 0
        }else{
          ans$RT_class[curr_idx[length(curr_idx)]] <- "P"
        }

        #merge similar elements into an RT mean
        for(j in 2:length(curr_idx)){
          #if current is P and previous is N and their difference is within the limit
          if(ans$RT_class[curr_idx[j]] == "P" & ans$RT_class[curr_idx[j-1]] == "N" & ((ans$RT[curr_idx[j]] - ans$RT[curr_idx[j-1]]) < rt_class) ){
            #extend curr_idx[j] to include all elements that have same RT
            temp <- tail(which(ans$Peak.Name_RT_mean == ans$Peak.Name_RT_mean[curr_idx[j]]), 1)
            #from previous to all elments with same RT as the current will get the mean retention time
            ans$RT_mean[curr_idx[j-1]:temp] <- mean(ans$RT_mean[curr_idx[j-1]:temp]) #update all elments between j-1 and j: 2, 5 -> 2, 3, 4, 5, 6
          }
        }

      }

    }

    #print(ans)

    #eliminate duplicated due to RT aggregation
    #ans <- unique(ans[, c("Peak.Name", "RT_mean")])
    n_classes_curr <- nrow(unique(ans[, c("Peak.Name", "RT_mean")]))

  }

  #ans <- merge(X, ans, by=c("Peak.Name"))
  ans <- ans[order(ans$Peak.Name, ans$RT_mean, ans$RT), c("Peak.Name", "RT", "Peak.Name_RT_mean", "RT_mean")]

  return(ans)

}
