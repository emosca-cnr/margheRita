#' PCA
#' Function that allows to perform PCA. Graphs of Scree plot and pairs of the first 10 components were plotted. Score and loading plot, choosing pcx and pcy, were plotted launching function plot2DPCA after pca function.
#' @param mRList mRList object
#' @param col_by (default class)
#' @param scaling choose the scaling method "none","Pareto" or "uv"
#' @param include_QC (default TRUE)
#' @param write_output (default =FALSE) if it turns on TRUE tables as .csv of score and loading will be saved
#' @return  mRList object with "pca" element
#' @export
#' @importFrom graphics plot barplot pairs
#' @importFrom grDevices dev.off png rainbow
#' @importFrom utils write.csv
#' @importFrom stats var prcomp
#' @param top only the top most varying features will be used
#' @param dirout output directory
#' @param rank number of PCs

pca_gen <- function(mRList, dirout, col_by="class", scaling=c("none", "Pareto", "uv"), include_QC=TRUE, top=Inf, write_output=FALSE, rank=10) {
  
  
  scaling <- match.arg(scaling)
  
  dir.create(dirout)
  
  X <- mRList$data
  X_ann <- mRList$sample_ann
  
  #include QC
  if(include_QC){
    cat("Including QC\n")
    X <- cbind(X, mRList$QC)
    X_ann<- rbind(X_ann, mRList$QC_ann)
  }
  
  if(nrow(X) > top){
    cat("using only top", top, "metabolites by variance\n")
    idx <- order(-apply(X, 1, var))[1:top]
    X <- X[idx, ]
  }
  
  scale <- FALSE
  center <- FALSE
  
  #pareto scaling
  if (scaling=="Pareto") {
    cat("Pareto scaling\n")
    X <- pareto(X) #ettore: pareto requires mRList
  }
  if(scaling=="uv"){
    cat("UV scaling\n")
    scale <- TRUE
    center <- TRUE
  }
  
  
  #it is the same as PCA for euclidean distance
  
  pca <- prcomp(X, scale = scale, center = center, rank=rank)
  p.v. <- matrix(((pca$sdev ^ 2) / (sum(pca$sdev ^ 2))), ncol = 1) #varianza
  p.i. <- round(p.v.* 100, 1) #percentuali di varianza spiegata dalle PC
  Pvar. <- p.i.
  #mRList$pca <- append(pca, list (variance=Pvar.))
  
  if (write_output){
    pwd.score= paste(dirout, "/ScoreMatrix.csv", sep ="")
    utils::write.csv(pca$x, pwd.score)
    pwd.load. = paste(dirout ,"/LoadingsMatrix.csv", sep= "")
    utils::write.csv(pca$rotation, pwd.load.)
    pwd.pvar.= paste(dirout, "/Variance.csv", sep = "")
    utils::write.csv(p.i., pwd.pvar.)
    
  }
  
  
  #Graphic pairs
  
  pairplot= paste(dirout, "/First_10_components.png", sep = "")
  grDevices::png(
    pairplot,
    width = 8,
    height = 8,
    units = "in",
    res= 300
  )
  
  col_factor <- as.factor(X_ann[, col_by])
  col_pal <- rainbow(length(levels(col_factor)))
  pairs(mRList$pca$x[, 1:10], labels=  Pvar.,
        main = "First 10 components",
        col = col_pal[as.numeric(col_factor)],
        pch = 19
  )
  
  grDevices::dev.off()
  
  
  #Graphic screeplot
  
  
  scree = paste(dirout, "/Screeplot.png", sep = "")
  grDevices::png(
    scree,
    width = 8,
    height = 8,
    units = "in",
    res = 300
  )
  graphics::barplot(
    Pvar.[1:rank, 1],
    xlab = "Principal Components",
    ylab = "Proportion of Variance explained",
    main = "Screeplot",
    ylim = c(0, 100)
  )
  grDevices::dev.off()
  
  
  return(mRList)
  
}
