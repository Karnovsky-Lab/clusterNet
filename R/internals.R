#' getConsensusMatrix will calculate the consensus matrix of the network
#'
#' The function takes as input the results from the consensus clustering algorithm and calculates
#' the consensus matrix
#'
#' @param cl the results from the consensus clustering algorithm
#' @return The corresponding consensus matrix
#' @keywords internal
getConsensusMatrix <- function(cl){
  K <- length(cl)
  p <- length(cl[[1]]$membership)
  D <- matrix(0, p, p)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      for (k in 1:K){
        tmp = cl[[k]]$membership
        D[i,j] <- D[i,j] + (length(unique(tmp[c(i,j)]))==1)
      }
    }
  }
  D <- D + t(D)
  D <- D/K
  diag(D) <- rep(1, p)
  return(D)
}
# catch_cluster_leading_eigen <- function(graph, weights = NULL){
#
#   clustering_results <-trycatch(
#     cluster_leading_eigen(graph, weights = weights),
#     error = function(e){
#       message('cluster_leading_eigen() method failed and will be discarded from consensus clustering.')
#       message('This is a known issue with a dependency and will not affect your results')
#       return(NA)
#     }
#   )
# }

#' Performs consensus clustering
#'
#' This function will take as input an adjacency matrix graph from the determined networks and perform
#' consensus clustering using the following methods from the igraph package: cluster_edge_betweenness,
#' cluster_fast_greedy, cluster_infomap, cluster_label_prop, cluster_leading_eigen, cluster_louvain,
#' cluster_walktrap. The output results in sub-network classification for the nodes within the network.
#'
#' @param graph An adjacency matrix of the determined network
#' @param K The number of iterations for clustering. This parameter is not necessary for default "ensemble"
#' @param tau The consensus probabilty threshold for agreement among clustering runs
#' @param method The method to use for consensus cluster
#' @param maxIter Maximum number of iterations to perform
#'
#' @return Sub-network determinations for the nodes within the input network
#'
#' @import igraph
#' @keywords internal
run_consensus_cluster <- function(adjacency_graph,
                                  num_iterations=10,
                                  tau=0.5,
                                  maxIter=5){

  clustering_results <- list()
  set.seed(1)

  clustering_results[[1]] <- cluster_edge_betweenness(adjacency_graph, weights = NULL)
  clustering_results[[2]] <- cluster_fast_greedy(adjacency_graph, weights = NULL)
  clustering_results[[3]] <- cluster_infomap(adjacency_graph, e.weights = NULL)
  clustering_results[[4]] <- cluster_label_prop(adjacency_graph, weights = NULL)
  clustering_results[[5]] <- cluster_louvain(adjacency_graph, weights = NULL)
  clustering_results[[6]] <- cluster_walktrap(adjacency_graph, weights = NULL)
  clustering_results[[7]] <- tryCatch(cluster_leading_eigen(adjacency_graph, weights = NULL),
                                      error = function(some_error){
                                        message('cluster_leading_eigen() method failed and will be discarded from consensus clustering.')
                                        message('This is a known issue with a dependency and will not affect your results')
                                        return(NA)
                                      })

  D <- getConsensusMatrix(clustering_results[!(is.na(clustering_results))])

  iter <- 0
  while(length(table(D))>1 && iter<maxIter){
    diag(D) <- 0
    thresholded_D <- D*(D>tau)

    Dgraph <- graph.adjacency(thresholded_D, mode="undirected", weighted = TRUE)

    dcl <- list()
    set.seed(iter)
    dcl[[1]] <- cluster_edge_betweenness(Dgraph, weights = E(Dgraph)$weight)
    dcl[[2]] <- cluster_fast_greedy(Dgraph, weights = E(Dgraph)$weight)
    dcl[[3]] <- cluster_infomap(Dgraph, e.weights = E(Dgraph)$weight)
    dcl[[4]] <- cluster_label_prop(Dgraph, weights = E(Dgraph)$weight)
    dcl[[5]] <- cluster_louvain(Dgraph, weights = E(Dgraph)$weight)
    dcl[[6]] <- cluster_walktrap(Dgraph, weights = E(Dgraph)$weight)
    dcl[[7]] <- tryCatch(cluster_leading_eigen(Dgraph,
                                               weights = E(Dgraph)$weight),
                         error = function(some_error){
                           message('cluster_leading_eigen() method failed and will be discarded from consensus clustering.')
                           message('This is a known issue with the igraphf package and will not affect your results')
                           return(NA)
                         })
    D <- getConsensusMatrix(dcl[!(is.na(dcl))])
    iter <- iter + 1
  }

  new.order <- order(dcl[[1]]$membership)
  return(list(dcl=dcl[[1]]$membership,D=D,order=new.order,iter=iter))
}
