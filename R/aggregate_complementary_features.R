aggregate_complementary_features <- function(X, rt_class=NULL){
  
  #"Peak.Name_RT_mean"	...
  n_samples <- ncol(X) - 1
  peaks_rt <- as.data.frame(do.call(rbind, strsplit(X[, 1], "_")), stringsAsFactors = F) #name, RT_mean
  peaks_rt[, 2] <- as.numeric(peaks_rt[, 2])
  
  #calculate the complementarity structure
  ans <- split(X, peaks_rt[, 1])
  peaks_rt <- split(peaks_rt[, 2], peaks_rt[, 1])
  
  table_RT <- vector("list", length(ans))
  names(table_RT) <- names(ans)
  table_compl <- table_RT
  table_compl_flag <- table_RT
  table_RT_flag <- table_RT
  for(i in 1:length(ans)){
    n_elements <- nrow(ans[[i]])
    
    table_RT[[i]] <- matrix(0, nrow = n_elements, ncol = n_elements)
    table_RT_flag[[i]] <- 0
    table_compl[[i]] <- matrix(0, nrow = n_elements, ncol = n_elements)
    table_compl_flag[[i]] <- 0
    if(n_elements>1){
      for(j in 2:n_elements){
        for(k in 1:(j-1)){
          table_RT[[i]][j, k] <- abs(peaks_rt[[i]][j] - peaks_rt[[i]][k])
          if(table_RT[[i]][j, k] < rt_class){
            table_RT_flag[[i]] <- table_RT_flag[[i]] + 1
          }
          table_compl[[i]][j, k] <- sum(colSums(sign(ans[[i]][c(j, k), -1])) == 1) #number of 1's
          complementary_pairs <- sum(colSums(sign(ans[[i]][c(j, k), -1])) < 2) / n_samples #number of 1's and 0's
          if(table_compl[[i]][j, k] > 0 & complementary_pairs == 1){ #not all 0's
            table_compl_flag[[i]] <- table_compl_flag[[i]] + 1
          }
        }
      }
    }
  }
  
 
  #aggregable
  #cat("aggregables...")
  #print(table(unlist(table_RT_flag) + unlist(table_compl_flag)))
  
  idx_aggr <- which(sign(unlist(table_RT_flag)) + sign(unlist(table_compl_flag)) == 2)
  
  if(length(idx_aggr) > 0){
    
    cat("Aggregation of:")
    print(names(idx_aggr))
    
    for(i in 1:length(idx_aggr)){
      new_row <- ans[[idx_aggr[i]]]
      
      idx_pair <- which(table_RT[[idx_aggr[i]]] < rt_class & table_RT[[idx_aggr[i]]] > 0, arr.ind = T) # indices of the aggregable pair(s)
      
      if(nrow(idx_pair) > 1){
        
        warning("multiple choices in aggregation of RT\n")
        cat("multiple choices in aggregation of RT\n")
        print(new_row[, 1:5])
        print(table_RT[[idx_aggr[i]]])
        print(table_compl[[idx_aggr[i]]])
        
        #selection of the one with the highest complementarity
        idx_pair <- idx_pair[which.max(table_compl[[idx_aggr[i]]][idx_pair]), ]
        cat("selected: ", idx_pair, "\n")

      }
    
      new_row <- t(apply(new_row[idx_pair, -1], 2, max)) #aggregation
      new_row <- data.frame(Peak.Name_RT=paste0(names(peaks_rt)[idx_aggr[i]], "_", mean(peaks_rt[[idx_aggr[i]]][idx_pair])), new_row, stringsAsFactors = F)
      rownames(new_row) <- new_row$Peak.Name_RT_mean
      ans[[idx_aggr[i]]] <- ans[[idx_aggr[i]]][-idx_pair, ]
      ans[[idx_aggr[i]]] <- rbind(ans[[idx_aggr[i]]], new_row)
    }
  }
  
  ans <- do.call(rbind, ans)
  rownames(ans) <- ans$Peak.Name_RT
  
  return(ans)
  
}
