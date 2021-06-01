#' Get connexted comnponents of size at least n
#' @import igraph
#' @export
#'
get.n.conn.comp <- function(x, n){

  conn <- igraph::clusters(x)
  clstr.ok <- which(conn$csize >= n)
  nodes.ok <- which(conn$membership %in% clstr.ok)

  return(igraph::induced.subgraph(x, nodes.ok))
}
