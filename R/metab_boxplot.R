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

  data<-m_list$data[select=features,] #select the metabolites (features), according to list from users
  data<-t(data)
  groups<-as.factor(X_ann[,group])


    for (i in 1:ncol(data)) {
      boxplot(
     data[,i] ~ groups , data=data,
      names =levels(groups),
      #main = m_list$metabo_ann, #metabolites or Ms ID or annotation
      ylab="Relative Abundance", col= col_pal)
    png(file=paste(dirout,"metabolite",i,".png",sep=""), width = 480, height = 480, units = "px",
        bg = "white")
  }
  dev.off()
 }
