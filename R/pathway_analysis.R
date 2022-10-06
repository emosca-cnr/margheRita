#' Pathway analysis
#' Run enricher or gsea from clusterProfiler on pathways provided by NCBI Biosytems
#' @param in_list metabolite set
#' @param universe universe
#' @param type "ora" or "gsea"
#' @param tax_id tax id
#' @param include_general_pathways whether to include general pathways or not
#' @param pvalueCutoff p-value cutoff
#' @param pAdjustMethod p-value adjustment method
#' @param minGSSize minimum metabolite set size
#' @param maxGSSize maximum metabolite set size
#' @param qvalueCutoff q-value cutoff
#' @param TERM2METAB 2-columns data.frame of pathways and metabolites
#' @param TERM2NAME 2-columns data.frame of pathways and pathway names
#' @param nPerm number of permutations
#' @param verbose verbose mode
#' @param gsea_by gsea type
#' @param seed set the seed
#' @param exponent exponent in the gsea approach
#' @export
#' @importFrom  clusterProfiler GSEA enricher
#' @importFrom utils data

pathway_analysis <- function(in_list=NULL, universe=NULL, type=c("ora", "gsea"), tax_id=9606, include_general_pathways=FALSE, pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, TERM2METAB=NULL, TERM2NAME = NULL, nPerm=1000, verbose=TRUE, gsea_by="fgsea", seed=FALSE, exponent=1){

  bsid2cid <- bsid2info <- NULL #to please the check
  
  type <- match.arg(type, c("ora", "gsea"))

  if(is.null(TERM2METAB)){

    cat("loading pathway data...\n")

    data("bsid2cid", envir=environment())
    data("bsid2info", envir=environment())

    cat("preprocessing pathway data...\n")

    #pathway definitions only from the given tax
    if(include_general_pathways){
      TERM2NAME <- bsid2info[which(bsid2info$tax_id == tax_id | is.na(bsid2info$tax_id)), ]
    }else{
      TERM2NAME <- bsid2info[which(bsid2info$tax_id == tax_id), ]
    }

    #pathway to metab only from the considered pathways
    TERM2METAB <- bsid2cid[bsid2cid$bs_id %in% TERM2NAME$bs_id, ]

    #only required columns
    TERM2METAB <- TERM2METAB[, c("bs_id", "CID")]
    TERM2NAME_ <- TERM2NAME #full annotation
    TERM2NAME <- TERM2NAME[, c("bs_id", "name")]

  }

  #restrict the pathway db to the given universe
  if(!is.null(universe)){

    cat("reducing pathway data to the given universe...\n")

    TERM2METAB <- TERM2METAB[TERM2METAB$CID %in% universe, ]

    TERM2NAME <- TERM2NAME[TERM2NAME$bs_id %in% TERM2METAB$bs_id, ]

    TERM2NAME_ <- TERM2NAME_[TERM2NAME_$bs_id %in% TERM2METAB$bs_id, ]
  }


  if(type=="ora"){

    cat("ORA...\n")
    res <- clusterProfiler::enricher(in_list, pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, universe=universe, minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff, TERM2GENE=TERM2METAB, TERM2NAME = TERM2NAME)


  }else{

    cat("GSEA...\n")

    res <- clusterProfiler::GSEA(geneList = in_list, exponent = exponent, nPerm = nPerm, minGSSize = minGSSize, maxGSSize = maxGSSize, pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, TERM2GENE = TERM2METAB, TERM2NAME = TERM2NAME, verbose = verbose, seed = seed, by = gsea_by)

  }

  res@result <- res@result
  #colnames(res@result) <- gsub("[Gg]ene", "Compound", colnames(res@result))

  return(list(res=res, term_annotation=TERM2NAME_))

}
