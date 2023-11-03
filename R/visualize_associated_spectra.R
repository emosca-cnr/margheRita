#' MS/MS Spectra Visualization
#' Visualize the MS/MS spectra of the features associated with a metabolite.
#'
#' @param mRList mRList object that contains the results of metabolite identification
#' @param mR_library library list used during annotation
#' @param metabolite_id ID of the metabolite to be visualized
#' @param out_dir output directory
#' @param peak_number if TRUE write the peak number
#' @export
#' @importFrom graphics par lines legend abline
#' @importFrom grDevices jpeg adjustcolor
#' @importFrom plotrix thigmophobe.labels
#' @importFrom Hmisc minor.tick
#' @param type mirrored or overlapped visualization

visualize_associated_spectra <- function(mRList=NULL, mR_library=NULL, metabolite_id=NULL, out_dir="./", type=c("mirrored", "overlapped"), peak_number=FALSE){

  type <- match.arg(type)

  metabolite_name <- unique(mRList$metabolite_identification$associations$Name[mRList$metabolite_identification$associations$ID == metabolite_id])
  
  if(length(metabolite_name) == 0){
    stop("Can't find ", metabolite_id, " between level 1-2 associations.\n")
  }
  if(!any(mRList$metabolite_identification$associations$Level[mRList$metabolite_identification$associations$Name %in% metabolite_name] %in% c(1, 2))){
    stop("Can't find ", metabolite_id, " between level 1-2 associations.\n")
  }

  #get the features associated with metabolite_id
  associated_features <- unique(mRList$metabolite_identification$associations[mRList$metabolite_identification$associations$ID == metabolite_id & mRList$metabolite_identification$associations$Level %in% c(1, 2), c("Feature_ID", "ID_peaks"), drop=F])

  RI_sample <- mRList$metabolite_identification$RI_sample
  RI_lib <- mR_library$lib_peaks

  for(i in 1:nrow(associated_features)){

    feature_id <- associated_features$Feature_ID[i]
    peak_id <- associated_features$ID_peaks[i]

    feature_spectra <- RI_sample[[as.character(feature_id)]]
    library_spectra <- RI_lib[[as.character(peak_id)]]
    match_matrix <- mRList$metabolite_identification$MS_MS_info[[metabolite_id]][[peak_id]][[as.character(feature_id)]]$flags

    lib_n <- which(apply(match_matrix, 1, function(x) any(x==3)))

    dir.create(out_dir, recursive = T)

    jpeg(filename = paste0(out_dir, "/", metabolite_id, "_", feature_id, "_",  peak_id, "_", type, ".jpg"), width = 200, height = 200, res=300, units = "mm")

    if(type=="mirrored"){
      ylim <- c(-100, 100)
      hlines <- seq(-100, 100, 10)
    }else{
      ylim <- c(0, 100)
      hlines <- seq(0, 100, 10)
    }

    #empty plot
    plot(0, pch="", xlim=c(min(library_spectra[, 1], feature_spectra[, 1]), max(library_spectra[, 1], feature_spectra[, 1])), ylim = ylim, xlab = "m/z", ylab = "Relative Intensity")
    abline(h=hlines, lty=2, col="gray")
    minor.tick(ny = 2, nx=1)

    #library
    points(library_spectra[, 1] , library_spectra[, 2], type = "h" , col= adjustcolor("red", 0.6), lwd=2)
    points(library_spectra[, 1], library_spectra[, 2], col=adjustcolor("red", 0.8), pch=16)

    if(peak_number){
      thigmophobe.labels(library_spectra[, 1], library_spectra[,2], 1:nrow(library_spectra))
    }

    ### feature
    if(type=="mirrored"){

      lines(feature_spectra[, 1], -feature_spectra[,2], type = "h" , col=adjustcolor("blue", 0.6), lwd=2)
      points(feature_spectra[, 1], -feature_spectra[,2], col=adjustcolor("blue", 0.6), pch=16)
      if(peak_number){
        thigmophobe.labels(feature_spectra[, 1], -feature_spectra[, 2], 1:nrow(feature_spectra))
      }

    }else{

      lines(feature_spectra[, 1], feature_spectra[,2], type = "h" , col=adjustcolor("blue", 0.6), lwd=2)
      points(feature_spectra[, 1], feature_spectra[,2], col=adjustcolor("blue", 0.6), pch=16)
      if(peak_number){
        thigmophobe.labels(feature_spectra[, 1], feature_spectra[, 2], 1:nrow(feature_spectra))
      }

    }


    if(!is.null(match_matrix)){
      lib_n <- which(apply(match_matrix, 1, function(x) any(x==3)))
      points(library_spectra[lib_n, 1], rep(0, length(lib_n)), pch=1, col="purple", cex=2)
    }

    legend("bottomright", legend=c(metabolite_name, feature_id, "match"), lty=c(1, 1, NA), col=c("red", "blue", "purple"), cex=0.7, pch=c(NA, NA, 1), bg="transparent")

    dev.off()

  }

}
