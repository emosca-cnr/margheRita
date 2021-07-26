#' margherRita - main function to run the full pipeline
#'
#'
#' @import notame MSnbase

margheRita_test_pa <- function(wdir="./"){

  input_data_file <- system.file("extdata", "dataset_drift.xlsx", package = "margheRita")
  input_metadata_file <- system.file("extdata", "dataset_drift_metadata.xlsx", package = "margheRita")

  ### 1 ### READ INPUT
  cat("reading:\n", input_data_file, "\n", input_metadata_file, "\n")
  m_list_init <- read_input_file(input_data_file, metadata = input_metadata_file, data_start_col = 8, rt_col = 2, mz_col = 3)
  lapply(m_list_init, head)


  ### 3 ### FILTER BY m/z
  m_list <- m_z_filtering(m_list_init) ### issue with m_z_average

  ### 4 ### FILTER BY MISSING VALUES
  m_list <- filter_NA(m_list)


  ### 9 ### NORMALIZATION
  m_list <- calc_reference(m_list)
  norm_data <- normalize_profiles(m_list, method = "pqn")
  rla_res <- RLA(norm_data, do_plot = T, outline=F, las=2, out_dir = "RLA_norm", pars=list(cex.axis=0.3))

  ## COLLASSO
  norm_data_biorep <- collapse_tech_rep(norm_data, remove.QC = FALSE)

  ### 7 ### FILTER BY CV
  norm_data_biorep <- CV(norm_data_biorep, dirout = "CV") #this is too slow

  # ## 11 ### CHECK for BATCH EFFECT


  #### PATHWAY ANALYSIS ####
  #toy metab list and ranked list to test pathway analysis
  bsid2cid <- bsid2cid[bsid2cid$CID %in% metabolite_annotation$PubChem.CID, ]
  bsid2pubchem <- bsid2pubchem[bsid2pubchem$bs_id %in% bsid2cid$bs_id, ]

  #ORA Example
  metab_list <- unique(c(bsid2cid$CID[bsid2cid$bs_id=="82926"], sample(unique(bsid2cid$CID), 10))) #glycolysis
  res_ora <- pathway_analysis(in_list = metab_list, universe = metabolite_annotation$PubChem.CID, minGSSize = 5)


  #GSEA Example
  metab_ranked_list <- unique(metabolite_annotation$PubChem.CID)
  metab_ranked_list <- setNames(length(metab_ranked_list):1, metab_ranked_list) #decreasing values
  max_val <- (max(metab_ranked_list)+100) #glycolysis genes at the top
  min_val <- (max(metab_ranked_list)+100-length(metab_list)+1) #glycolysis genes at the top
  metab_ranked_list[names(metab_ranked_list) %in% metab_list] <- max_val:min_val #glycolysis genes at the top
  metab_ranked_list <- sort(metab_ranked_list, decreasing = T) #glycolysis genes at the top
  res_gsea <- pathway_analysis(in_list = metab_ranked_list, type = "gsea", universe = metabolite_annotation$PubChem.CID, minGSSize = 5, nPerm=99, pvalueCutoff = 1)



}
