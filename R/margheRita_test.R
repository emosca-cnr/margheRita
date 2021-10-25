#' margherRita - main function to run the full pipeline
#'
#'
#' @import notame MSnbase

margheRita_test <- function(wdir="./"){

  input_data_file <- system.file("extdata", "example1.xlsx", package = "margheRita")
  input_metadata_file <- system.file("extdata", "example1_meta.xlsx", package = "margheRita")

  ### 1 ### READ INPUT
  cat("reading:\n", input_data_file, "\n", input_metadata_file, "\n")
  m_list_init <- read_input_file(input_data_file, metadata = input_metadata_file)
  lapply(m_list_init, head)

  #mset <- as.metaboset.mRList(mRlist)
  #diff_res <- notame::perform_pairwise_t_test(mset)

  #clustered <- cluster_features(mset$HILIC_pos, all_features = T)
  #compressed <- compress_clusters(clustered)

  ### create MSnSet
  #msnset <- as.MSnSet.mRList(mRlist)

  #read mgf library
  #gmf_file <- system.file("extdata", "PSU-MSMLS.mgf", package = "margheRita")
  #temp <- MSnbase::readMgfData(gmf_file)


  ### 2 ### PLOTS
  pca_gen(m_list_init, dirout = "pca_initial") #fix the group variable for coloring

  ### 3 ### FILTER BY m/z
  m_list <- m_z_filtering(m_list_init) ### issue with m_z_average

  ### 4 ### FILTER BY MISSING VALUES
  m_list <- filter_NA(m_list)

  ### 5 ### PLOTS
  pca_gen(m_list, dirout = "pca_filtered") #fix the group variable for coloring

  ### 6 ### IMPUTATION
  m_list <- imputation(m_list) #this is too slow...


  ### 8 ### PLOTS
  pca_gen(m_list, dirout = "pca_imp") #fix the group variable for coloring
  rla_res <- RLA(m_list, include_QC=TRUE, do_plot = T, outline=F, las=2, out_dir = "RLA_raw", pars=list(cex.axis=0.3))

  ### 9 ### NORMALIZATION
  m_list <- calc_reference(m_list, sample_col = )
  norm_data <- normalize_profiles(m_list, method = "pqn")
  rla_res <- RLA(norm_data, do_plot = T, outline=F, las=2, out_dir = "RLA_norm", pars=list(cex.axis=0.3))

  ### 10 ### PLOTS
  pca_gen(norm_data, dirout = "pca_norm") #fix the group variable for coloring
  pca_gen(norm_data, dirout = "pca_norm_batch", col_by = "batch") #fix the group variable for coloring

  ## COLLASSO
  norm_data_biorep <- collapse_tech_rep(norm_data, remove.QC = FALSE)
  norm_data_biorep_ <- mean_media_stdev_samples(norm_data_biorep, dirout = "")

  h_map(norm_data_biorep)

  ### 7 ### FILTER BY CV
  norm_data_biorep <- CV(norm_data_biorep, dirout = "CV") #this is too slow

  # ## 11 ### CHECK for BATCH EFFECT

  ### 12 ### BATCH EFFECT removal
  temp <- batch_effect(norm_data)
  pca_gen(temp, dirout = "pca_batch") #fix the group variable for coloring
  RLA(temp, do_plot = T, outline=F, las=2, out_dir = "RLA_batch", pars=list(cex.axis=0.3))

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
