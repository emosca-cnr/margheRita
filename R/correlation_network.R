#' correlation network
#' @export
correlation_network <- function(m_list, top=0.05, query_metab=NULL){

  #under construction...

  #compute the correlation matrix
  if(!is.null(query_metab)){

    idx_qry <- rownames(m_list$data) %in% query_metab
    res_cor <- cor(t(m_list$data[idx_qry, ]), t(m_list$data))

    #rows
    res_cor <- rbind(res_cor, matrix(0, nrow = ncol(res_cor)-length(query_metab), ncol = ncol(res_cor), dimnames = list(colnames(res_cor)[!idx_qry], colnames(res_cor))))

    res_cor <- res_cor[match(colnames(res_cor), rownames(res_cor)), ]
    isSymmetric(res_cor) #no because there are correlations between query and the others... (I think)
  }else{

    res_cor <- cor(t(m_list$data))

  }

  #select the most interesting links
  diag(res_cor) <- 0
  topq <- quantile(abs(res_cor), probs=1-top)

  res_cor_filtered <- res_cor
  res_cor_filtered[abs(res_cor_filtered) < topq] <- 0

  #define the igraph object
  g <- igraph::graph.adjacency(res_cor_filtered, mode = "undirected", weighted = T)
  print(g)

  #remove single nodes
  g <- get.n.conn.comp(g, 2)
  print(g)

  #calculate the communities
  g_comm <- igraph::fastgreedy.community(g)

  #calculate the centralities



}
