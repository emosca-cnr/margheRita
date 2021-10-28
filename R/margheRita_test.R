#' margherRita - main function to run the full pipeline
#'
#'
#' @import notame MSnbase

margheRita_test <- function(wdir="./"){

  input_data_file <- system.file("extdata", "example1.xlsx", package = "margheRita")
  input_metadata_file <- system.file("extdata", "example1_meta.xlsx", package = "margheRita")

  ### 1 ### READ INPUT
  cat("reading:\n", input_data_file, "\n", input_metadata_file, "\n")
  mRList_raw <- read_input_file(input_data_file, metadata = input_metadata_file)
  lapply(mRList_raw, head)

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
  pca_gen(mRList_raw, dirout = "pca_initial") #fix the group variable for coloring

  res <- metab_boxplot(mRList_raw, features = c("M13167", "M25188"))

  ### 3 ### FILTER BY m/z
  mRList_filt <- m_z_filtering(mRList_raw) ### issue with m_z_average

  ### 4 ### FILTER BY MISSING VALUES
  mRList_filt <- filter_NA(mRList_filt)

  ### 5 ### PLOTS
  pca_gen(mRList_filt, dirout = "pca_filtered") #fix the group variable for coloring

  ### 6 ### IMPUTATION
  mRList_filt <- imputation(mRList_filt) #this is too slow...


  ### 8 ### PLOTS
  pca_gen(mRList_filt, dirout = "pca_imp") #fix the group variable for coloring
  rla_res <- RLA(mRList_filt, include_QC=TRUE, do_plot = T, outline=F, las=2, out_dir = "RLA_raw", pars=list(cex.axis=0.3))

  ### 9 ### NORMALIZATION
  mRList_filt <- calc_reference(mRList_filt)
  mRList_norm <- normalize_profiles(mRList_filt, method = "pqn")
  rla_res <- RLA(mRList_norm, do_plot = T, outline=F, las=2, out_dir = "RLA_norm", pars=list(cex.axis=0.3))

  ### 10 ### PLOTS
  pca_gen(mRList_norm, dirout = "pca_norm") #fix the group variable for coloring
  pca_gen(mRList_norm, dirout = "pca_norm_batch", col_by = "batch") #fix the group variable for coloring

  ## COLLASSO
  mRList_norm_biorep <- collapse_tech_rep(mRList_norm, remove.QC = FALSE)
  mRList_norm_biorep_ <- mean_median_stdev_samples(mRList_norm_biorep, dirout = "")

  h_map(mRList_norm_biorep)

  ### 7 ### FILTER BY CV
  mRList_norm_biorep <- CV(mRList_norm_biorep, dirout = "CV")

  # ## 11 ### CHECK for BATCH EFFECT

  ### 12 ### BATCH EFFECT removal
  temp <- batch_effect(mRList_norm)
  pca_gen(temp, dirout = "pca_batch") #fix the group variable for coloring
  RLA(temp, do_plot = T, outline=F, las=2, out_dir = "RLA_batch", pars=list(cex.axis=0.3))


  ### UNIVARIATE

  mRList_norm_biorep_ <- univariate(mRList_norm_biorep)

  ### ANNOTATION
  data4annot <- generate_dataset_for_annotation()

  rt_res <- check_RT(data4annot$sample_data, data4annot$lib_data)
  ppm_res <- check_mass(data4annot$sample_data, data4annot$lib_data)

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
