#' metab_boxplot
#' @importFrom graphics boxplot
#' @param m_list m_list
#' @export
#' @importFrom grDevices dev.off png

metab_boxplot<-function(m_list, dirout="./",features=NULL, col_by="class",group="class"){

  dirout = paste(dirout, sep = "")
  dir.create(dirout)

  X_ann <- m_list$sample_ann
  col_factor <- as.factor(X_ann[, col_by])
  col_pal <- rainbow(length(levels(col_factor)))

  data <- m_list$data[rownames(m_list$data) %in% features, ] #select the metabolites (features), according to list from users
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
