#' margherRita - main function to run the full pipeline
#'
#'
#' @import notame MSnbase

margheRita_test_pa <- function(wdir="./"){

  input_data_file <- system.file("extdata", "example1.xlsx", package = "margheRita")
  input_metadata_file <- system.file("extdata", "example1_meta.xlsx", package = "margheRita")

  ### 1 ### READ INPUT
  cat("reading:\n", input_data_file, "\n", input_metadata_file, "\n")
  m_list_init <- read_input_file(input_data_file, metadata = input_metadata_file, data_start_col = 4, rt_col = 2, mz_col = 3)
  lapply(m_list_init, head)

  ### 3 ### FILTER BY m/z
  m_list <- m_z_filtering(m_list_init) ### issue with m_z_average

  ### 4 ### FILTER BY MISSING VALUES
  m_list <- filter_NA(m_list)

  ### 9 ### NORMALIZATION
  m_list <- calc_reference(m_list)
  norm_data <- normalize_profiles(m_list, method = "pqn")
  rla_res <- RLA(norm_data, do_plot = T, outline=F, las=2, out_dir = "../../PROMEFA/margheRita/test/RLA_norm", pars=list(cex.axis=0.3))

  ## COLLASSO
  norm_data_biorep <- collapse_tech_rep(norm_data, remove.QC = FALSE)

  ### 7 ### FILTER BY CV
  norm_data_biorep <- CV(norm_data_biorep, dirout = "CV") #this is too slow

  rla_res <- RLA(norm_data_biorep, do_plot = T, outline=F, las=2, out_dir = "../../PROMEFA/margheRita/test/RLA_norm_collapsed", pars=list(cex.axis=0.3))


  ttest_res <- apply(norm_data_biorep$data[, grepl("STAND", colnames(norm_data_biorep$data))], 1, function(x) t.test(x ~ norm_data_biorep$sample_ann$class[grepl("STAND", norm_data_biorep$sample_ann$class)]))

  ttest_res_df <- data.frame(Feature_ID=names(ttest_res), t=unlist(lapply(ttest_res, function(x) x$statistic)), p=unlist(lapply(ttest_res, function(x) x$p.value)), stringsAsFactors = F)
  ttest_res_df$q <- p.adjust(ttest_res_df$p)

  mlist_sig <- data.frame(Feature_ID=ttest_res_df$Feature_ID[ttest_res_df$p<0.0005], stringsAsFactors = F)
  mlist_sig <- cbind(mlist_ann, CID=sample(metabolite_annotation$PubChem.CID, nrow(mlist_ann)))

  #### PATHWAY ANALYSIS ####

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
