#missing values imputation
#@param dataframe of input contains intensity of all features according to samples

imputation<- function(df) {
  for (i in 1:nrow(df[, 4:ncol(df)])) {
    set.seed(1234)
    min = min(df[i,],na.rm=T)
    data[i, is.na(df[i, ])] <-
      sample(seq(
        from = (min * 0.1),
        to = (min * 0.25),
        by = (min * 0.01)
      ), size = sum(is.na(df[i, ])), replace = T)
  }
  return(df)
}

