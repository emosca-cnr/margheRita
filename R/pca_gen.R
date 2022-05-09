#' PCA general
#' Scree plot, pairs of the first 10 components
#' Score and loading plot choosing pcx and pcy with function Plot2DPCA after pca_gen
#' @param mRList margheRita mRList
#' @param include_QC (default TRUE); scaling: Pareto, none,UV;col_by (default class)
#' @param write_output (default =FALSE) if it turns on TRUE tables as .csv of score and loading will be saved
#' @return Graphs of screeplot, pairs, scoreplot and loading plot; table as .csv if write_output turns TRUE
#' @export
#' @importFrom graphics plot barplot
#' @importFrom grDevices dev.off png
#' @importFrom utils write.csv


        #PCA function general to use in different points

        pca_gen <- function(mRList, dirout, col_by="class", scaling=c("none", "Pareto", "uv"),include_QC=TRUE, type=c("component", "coordinate"), dist.method="euclidean", top=Inf, write_output=FALSE) {


          type <- match.arg(type)
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

          if(type=="component"){

            #it is the same as PCA for euclidean distance

            pca <- prcomp(t(X), scale = scale, center = center)
            p.v.= matrix(((pca$sdev ^ 2) / (sum(pca$sdev ^ 2))), ncol = 1) #varianza
            p.i. = round(p.v.* 100, 1) #percentuali di varianza spiegata dalle PC
            Pvar. = p.i.
            mRList$pca<-append(pca, list (variance=Pvar.))
          }


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
       Pvar.[, 1],
       xlab = "Principal Components",
       ylab = "Proportion of Variance explained",
       main = "Screeplot",
       ylim = c(0, 100)
     )
     grDevices::dev.off()


           return(mRList)
           return(pairplot)
           return(scree)
           return(col_factor)
           return(col_pal)
        }



#'Create 2D PCA score plot
#'@param mRList Input name of the created object
#'@param scoreplot input a name for the score plot
#'@param loadingplot input for the loading plot
#'@param format Select the image format, "png".
#'@param pcx Specify the principal component on the x-axis
#'@param pcy Specify the principal component on the y-axis
#'@export
#'
Plot2DPCA <- function(mRList, pcx, pcy,dirout,col_by="class",include_QC=T){

  dir.create(dirout)

          xlabel = paste("PC",pcx, "(", mRList$pca$variance[pcx],1, "%)")
          ylabel = paste("PC",pcy, "(", mRList$pca$variance[pcy],1, "%)")
          pc1 = mRList$pca$x[, pcx]
          pc2 = mRList$pca$x[, pcy]

          X <- mRList$data
          X_ann <- mRList$sample_ann

          #include QC
          if(include_QC){
            cat("Including QC\n")
            X <- cbind(X, mRList$QC)
            X_ann<- rbind(X_ann, mRList$QC_ann)
          }


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
            pc1,pc2,
            xlab=xlabel,ylab=ylabel,
            main = "Score plot",
            col = col_pal[as.numeric(col_factor)],
            pch = 19
          )
          legend("bottomright", legend = levels(col_factor), col = col_pal, pch=16, cex=0.5)

          grDevices::dev.off()


          xlabeL = paste("Loading",pcx)
          ylabeL = paste("Loading",pcy)
          pc1L = mRList$pca$rotation[, pcx]
          pc2L = mRList$pca$rotation[, pcy]

          loadingplot = paste(dirout, "/Loadingplot.png", sep = "")
          grDevices::png(
            loadingplot,
            width = 8,
            height = 8,
            units = "in",
            res = 300
          )
          graphics::plot(
            pc1L,pc2L,
            xlab=xlabeL,ylab=ylabeL,
            main = "Loading plot",
            pch = 19
          )
          grDevices::dev.off()
}
