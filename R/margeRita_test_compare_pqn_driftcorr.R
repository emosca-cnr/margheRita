#' margherRita - main function to run the full pipeline
#'
#'
#' @import notame MSnbase

margheRita_test_compare_pqn_driftcorr <- function(wdir="./"){

  input_data_file <- system.file("extdata", "dataset_drift.xlsx", package = "margheRita")
  input_metadata_file <- system.file("extdata", "dataset_drift_metadata.xlsx", package = "margheRita")

  ### 1 ### READ INPUT
  cat("reading:\n", input_data_file, "\n", input_metadata_file, "\n")
  m_list_init <- read_input_file(input_data_file, metadata = input_metadata_file, data_start_col = 8, rt_col = 2, mz_col = 3)
  lapply(m_list_init, head)

  rla_res <- RLA(m_list_init, do_plot = T, outline=F, las=2, out_dir = "../../PROMEFA/margheRita/test/RLA_raw", pars=list(cex.axis=0.3))

  ### 3 ### FILTER BY m/z
  m_list <- m_z_filtering(m_list_init) ### issue with m_z_average

  ### 4 ### FILTER BY MISSING VALUES
  m_list <- filter_NA(m_list)

  ### 9 ### NORMALIZATION
  m_list <- calc_reference(m_list)
  norm_data <- normalize_profiles(m_list, method = "pqn")
  rla_res <- RLA(norm_data, do_plot = T, outline=F, las=2, out_dir = "../../PROMEFA/margheRita/test/RLA_norm", pars=list(cex.axis=0.3))

  ### 1.1 #### create a notame MetaboSet 0bject
  feature_data <- m_list_init$metab_ann
  colnames(feature_data)[3] <- "Average mz"
  feature_data <- cbind(feature_data, Split="HILIC_pos")
  feature_data$Feature_ID <- paste0("ID", feature_data$Feature_ID)
  rownames(feature_data) <- feature_data$Feature_ID

  temp <- cbind(m_list_init$data, m_list_init$QC)
  rownames(temp) <- feature_data$Feature_ID

  pheno_data <- rbind(m_list_init$sample_ann, m_list_init$QC_ann)
  colnames(pheno_data)[1] <- "Sample_ID"
  colnames(pheno_data)[2] <- "Injection_order"
  colnames(pheno_data)[4] <- "QC"

  mset <- notame::construct_metabosets(exprs = as.matrix(temp), pheno_data = pheno_data, feature_data = feature_data, group_col = "class")

  detected <- flag_detection(object, qc_limit = 0.7, group_limit = 0.5)
  mset <- notame::correct_drift(mset$HILIC_pos)

  mlist_notame_drift <- list(data = exprs(mset))
  rla_notame <- RLA(mlist_notame_drift, do_plot = T, outline=F, las=2, out_dir = "../../PROMEFA/margheRita/test/RLA_drift_notame", pars=list(cex.axis=0.3))


}
