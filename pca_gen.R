#' PCA general
#' Scree plot, score and loading plot (PC1 vs PC2)
#' @param m_list margheRita mRList
#' @param include_QC (default TRUE); scaling: Pareto, none,UV;col_by (default class)
#' @param write_output (default =FALSE) if it turns on TRUE tables as .csv of score and loading will be saved
#' @return Graphs of screeplot, scoreplot and loading plot; table as .csv if write_output turns TRUE
#' @export
#' @importFrom graphics plot barplot
#' @importFrom grDevices dev.off png
#' @importFrom utils write.csv


#PCA function general to use in different points

pca_gen <- function(m_list, dirout, col_by="class", scaling=c("none", "Pareto", "uv"), include_QC=TRUE, type=c("component", "coordinate"), dist.method="euclidean", top=Inf, write_output=FALSE) {


  type <- match.arg(type)
  scaling <- match.arg(scaling)

  dir.create(dirout)

  X <- m_list$data
  X_ann <- m_list$sample_ann

  #include QC
  if(include_QC){
    cat("Including QC\n")
    X <- cbind(X, m_list$QC)
    X_ann<- rbind(X_ann, m_list$QC_ann)
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
    X <- pareto(X) #ettore: pareto requires m_list
  }
  if(scaling=="uv"){
    cat("UV scaling\n")
    scale <- TRUE
    center <- TRUE
  }

  if(type=="component"){

    #it is the same as PCA for euclidean distance
    pca <- prcomp(t(X), scale = scale, center = center)

    p.v.= matrix(((pca$sdev ^ 2) / (sum(pca$sdev ^ 2))), ncol = 1) #varianza
    p.i. = round(p.v.* 100, 1) #percentuali di varianza spiegata dalle PC
    Pvar. = p.i.

    if (write_output){
    pwd.score= paste(dirout, "/ScoreMatrix.csv", sep ="")
    utils::write.csv(pca$x, pwd.score)
    pwd.load. = paste(dirout ,"/LoadingsMatrix.csv", sep= "")
    utils::write.csv(pca$rotation, pwd.load.)
    pwd.pvar.= paste(dirout, "/Variance.csv", sep = "")
    utils::write.csv(p.i., pwd.pvar.)

}

    #ora faccio i grafici di score, loading and scree plot plotting PC1 vs PC2

   scoreplot = paste(dirout, "/Scoreplot.png", sep = "")
    grDevices::png(
      scoreplot,
      width = 8,
      height = 8,
      units = "in",
      res = 300
    )

    col_factor <- as.factor(X_ann[, col_by])
    col_pal <- rainbow(length(levels(col_factor)))
    graphics::plot(
      pca$x[, 1],
      pca$x[, 2],
      xlab = paste("PC1 (", p.i.[1], "%)", sep = ""),
      ylab = paste("PC2 (", p.i.[2], "%)", sep = ""),
      main = "Score plot",
      col = col_pal[as.numeric(col_factor)],
      pch = 19
    )
    legend("bottomright", legend = levels(col_factor), col = col_pal, pch=16, cex=0.5)

    grDevices::dev.off()
    loadingplot = paste(dirout, "/Loadingplot.png", sep = "")
    grDevices::png(
      loadingplot,
      width = 8,
      height = 8,
      units = "in",
      res = 300
    )
    graphics::plot(
      pca$rotation[, 1],
      pca$rotation[, 2],
      xlab = "Loading 1",
      ylab = "Loading 2",
      main = "Loading plot",
      pch = 19
    )
    grDevices::dev.off()
    scree = paste(dirout, "/Screeplot.png", sep = "")
    grDevices::png(
      scree,
      width = 8,
      height = 8,
      units = "in",
      res = 300
    )
    graphics::barplot(
      Pvar.[, 1],
      xlab = "Principal Components",
      ylab = "Proportion of Variance explained",
      main = "Screeplot",
      ylim = c(0, 100)
    )
    grDevices::dev.off()

    }


  if(type=="coordinate"){

    pca <- cmdscale(dist(t(X), method = dist.method))

       scoreplot <- paste(dirout, "/Scoreplot.png", sep = "")
    grDevices::png(
      scoreplot,
      width = 8,
      height = 8,
      units = "in",
      res = 300
    )
    col_factor <- as.factor(X_ann[, col_by])
    col_pal <- rainbow(length(levels(col_factor)))
    graphics::plot(
      pca[, 1],
      pca[, 2],
      xlab = "PC1",
      ylab = "PC2",
      main = "Principal Coordinates",
      col = col_pal[as.numeric(col_factor)],
      pch = 19
    )
    legend("bottomright", legend = levels(col_factor), col = col_pal, pch=16, cex=0.5)
    dev.off()
  }

}
