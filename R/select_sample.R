#' sample selection
#'
#' @export
#' @param mRList obtained by read_input_file() function
#' @param n_column it is selecting the required column for metabolite annotation which are; Feature_ID, rt, mz, MS_MS_spectrum.
#'@param accept_RI numeric parameter. the default value is 10. it is a maximum relative intensity that is kept in sample. since low intense peaks could be noise, it is filtering sample data set by deleting the relative intensity lower then accept_RI.
#'
select_sample <- function(mRList = NULL, n_column = c(1, 2, 3, 5), accept_RI = 10) {
    sample <- mRList_raw$metab_ann[, n_column]
    sample_data <- as.data.frame(sample[!is.na(sample$MS_MS_spectrum), ])
    sample_data <- as.data.frame(sample_data[!is.na(sample_data$rt), ])
    sample_data <- as.data.frame(sample_data[!is.na(sample_data$mz), ])
    colnames(sample_data) <- c("Feature_ID", "rt", "mz", "MS_MS")


    sample_spectra <- sapply(sample_data$MS_MS, function(x) strsplit(x, " "))
    sample_spectra <- lapply(sample_spectra, function(x) lapply(strsplit(x, ":"), as.numeric))
    sample_spectra <- lapply(sample_spectra, function(x) do.call(rbind, x))
    names(sample_spectra) <- sample_data$Feature_ID
    sample_data$MS_MS <- NULL

    #converting intensity to relative intensity
    RI_sample <- RI_sample_data (feature_spectra = sample_spectra, accept_RI = accept_RI)

    return(list(feature_data = sample_data, feature_spectra = RI_sample))
  }
