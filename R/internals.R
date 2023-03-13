#' gatherConsensusMatrix will put together a model.matrix type summary of subnetwork results
#'
#' The function takes as input the metabolite info and subnetwork classfication to create an expanded matrix
#'
#' @param consensus_membership A numeric vector of subnetwork classification for each metabolite
#'
#' @returns a matrix with subnetworks as rows and metabolites as columns. The values are 0 if the metabolite is NOT
#'          in said subnetwork and a 1 if it is.
#' @keywords internal
gatherConsensusMatrix <- function(consensus_membership){

  subnetwork_results <- matrix(0, nrow=length(unique(consensus_membership)), length(consensus_membership))

  rownames(subnetwork_results) <- paste0("Subnetwork",1:length(unique(consensus_membership)))
  for (j in 1:nrow(subnetwork_results)){
    subnetwork_results[j,which(consensus_membership==j)] <- 1
  }
  if (length(which(rowSums(subnetwork_results)<5))>0){
    subnetwork_results <- subnetwork_results[-which(rowSums(subnetwork_results)<5),]
  }

  return(subnetwork_results)
}
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
#' @param adjacency_graph graph An adjacency matrix of the determined network.
#' @param tau The consensus probabilty threshold for agreement among clustering runs
#' @param num_iterations The total number of iterations to perform.
#' @param maxIter Maximum number of iterations to perform.
#' @param main.seed An integer to set seed for clustering. This parameter ensures reproduciblity and can be
#'                  kept as default.
#'
#' @return Sub-network determinations for the nodes within the input network
#'
#' @import igraph
#' @keywords internal
run_consensus_cluster <- function(adjacency_graph,
                                  tau,
                                  num_iterations=10,
                                  maxIter=5,
                                  main.seed){

  clustering_results <- list()
  if(!missing(main.seed)){

    set.seed(main.seed)
  } else{
    warning("Set seed to ensure reproducibility")
  }

  clustering_results[[1]] <- igraph::cluster_edge_betweenness(adjacency_graph, weights = NULL)
  clustering_results[[2]] <- igraph::cluster_fast_greedy(adjacency_graph, weights = NULL)
  clustering_results[[3]] <- igraph::cluster_infomap(adjacency_graph, e.weights = NULL)
  clustering_results[[4]] <- igraph::cluster_label_prop(adjacency_graph, weights = NULL)
  clustering_results[[5]] <- igraph::cluster_louvain(adjacency_graph, weights = NULL)
  clustering_results[[6]] <- igraph::cluster_walktrap(adjacency_graph, weights = NULL)
  clustering_results[[7]] <- tryCatch(igraph::cluster_leading_eigen(adjacency_graph, weights = NULL),
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
    if(!missing(main.seed)){

      set.seed(main.seed + iter)

    }else{

      set.seed(iter)

    }

    dcl[[1]] <- igraph::cluster_edge_betweenness(Dgraph, weights = E(Dgraph)$weight)
    dcl[[2]] <- igraph::cluster_fast_greedy(Dgraph, weights = E(Dgraph)$weight)
    dcl[[3]] <- igraph::cluster_infomap(Dgraph, e.weights = E(Dgraph)$weight)
    dcl[[4]] <- igraph::cluster_label_prop(Dgraph, weights = E(Dgraph)$weight)
    dcl[[5]] <- igraph::cluster_louvain(Dgraph, weights = E(Dgraph)$weight)
    dcl[[6]] <- igraph::cluster_walktrap(Dgraph, weights = E(Dgraph)$weight)
    dcl[[7]] <- tryCatch(igraph::cluster_leading_eigen(Dgraph,weights = E(Dgraph)$weight),
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
