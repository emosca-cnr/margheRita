# #Remove m/z > 0.4 or m/z < 0.8
# #@param dataframe input
# df$quality <- NA
# for(i in 1:nrow(df)) {
#   df$quality[i] <- (df$Average_mz[i] %% 1)
# }
# df <- filter(df, df$quality < 0.4 | df > 0.8 )
#
# plot(density(df$quality))
