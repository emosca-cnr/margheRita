#' Boxplot of metabolie levels across classes
#' 
#' Plot boxplot of metabolites, and save them as .png.
#'
#' @importFrom graphics boxplot
#' @param mRList mRList object
#' @param features vector with a list of metabolites to graph
#' @param col_by define how to color the boxplot, default=class
#' @param group define what you want to compare
#' @param dirout output directory
#' @export
#' @importFrom grDevices dev.off png

metab_boxplot<-function(mRList=NULL, dirout="./", features=NULL, col_by="class", group="class"){

  dirout = paste(dirout, sep = "")
  dir.create(dirout)

  X_ann <- mRList$sample_ann
  col_factor <- as.factor(X_ann[, col_by])
  col_pal <- rainbow(length(levels(col_factor)))

  data <- mRList$data[rownames(mRList$data) %in% features, ] #select the metabolites (features), according to list from users

  if(nrow(data)<1){
    message("None of the features was found in the dataset.\n")
  }

  #data<-t(data)
  groups<-as.factor(X_ann[,group])


  for (i in 1:nrow(data)) {

    png(file=paste0(dirout, "metabolite", rownames(data)[i], ".png"), width = 200, height = 200, units = "mm", res=300)

    i_min <- min(data[i, ])
    i_max <- max(data[i, ])

    graphics::boxplot(
      as.numeric(data[i, ]) ~ groups , data=data,
      names =levels(groups),
      main = rownames(data)[i], #metabolites or Ms ID or annotation
      ylab="Intensity",
      border= col_pal,
      ylim=c(i_min, i_max),
      outline=F,
      col="white"
    )
    for(j in 1:length(levels(groups))){
      idx_i <- which(groups==levels(groups)[j])
      points(jitter(rep(j, length(idx_i)), amount = 0.3), data[i, idx_i], col=col_pal[j], pch=16)
    }

    dev.off()

  }

}
