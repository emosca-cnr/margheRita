#' margherRita - main function to run the full pipeline
#'
#'

margheRita_test_pa <- function(wdir="./"){

 #### PATHWAY ANALYSIS ####
  bsid2cid <- bsid2cid[bsid2cid$CID %in% metabolite_annotation$PubChem.CID, ]
  bsid2info <- bsid2info[bsid2info$bs_id %in% bsid2cid$bs_id, ]

  #ORA Example
  metab_list <- unique(c(bsid2cid$CID[bsid2cid$bs_id=="82926"], sample(unique(bsid2cid$CID), 10))) #glycolysis
  res_ora <- pathway_analysis(in_list = metab_list, universe = metabolite_annotation$PubChem.CID, minGSSize = 5)

  jpeg(paste0(wdir, "/ora_barplot.jpg"), width = 180, height = 180, res=300, units="mm")
  barplot(res_ora$res, showCategory = 10, font.size = 12)
  dev.off()

  jpeg(paste0(wdir, "/ora_cnetplot.jpg"), width = 180, height = 180, res=300, units="mm")
  cnetplot(res_ora$res)
  dev.off()

  #GSEA Example
  metab_ranked_list <- unique(metabolite_annotation$PubChem.CID)
  metab_ranked_list <- setNames(length(metab_ranked_list):1, metab_ranked_list) #decreasing values

  #max_val <- (max(metab_ranked_list)+100) #glycolysis genes at the top
  #min_val <- (max(metab_ranked_list)+100-length(metab_list)+1) #glycolysis genes at the top

  max_val <- jitter(sample(metab_ranked_list[1:100], length(metab_list)), amount = 0.5) #glycolysis genes at the top

  metab_ranked_list[names(metab_ranked_list) %in% metab_list] <- max_val #glycolysis genes at the top
  metab_ranked_list <- sort(metab_ranked_list, decreasing = T) #glycolysis genes at the top

  res_gsea <- pathway_analysis(in_list = metab_ranked_list, type = "gsea", universe = metabolite_annotation$PubChem.CID, minGSSize = 5, nPerm=99, pvalueCutoff = 1)

  jpeg(paste0(wdir, "/gsea_runningscore.jpg"), width = 180, height = 180, res=300, units="mm")
  gseaplot(res_gsea$res, geneSetID = "82926", by = "runningScore")
  dev.off()

}
