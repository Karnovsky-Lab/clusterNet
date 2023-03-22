#' clusterNet performs consensus clustering using an adjacency matrix for the network
#'
#' This function will take as input an adjacency matrix from the determined network and perform
#' consensus clustering using the following methods from the igraph package: cluster_edge_betweenness,
#' cluster_fast_greedy, cluster_infomap, cluster_label_prop, cluster_leading_eigen, cluster_louvain,
#' cluster_walktrap. The output results in sub-network classification for the nodes within the network.
#'
#' @param data A network file, whether it be a correlation matrix or edge list
#' @param type a string identifying the type of input for the data parameter. **adjacency_matrix** is used for a
#'             correlation matrix and **edge_list** is used for inputting an edge list.
#' @param tau The consensus probabilty threshold for agreement among clustering algorithms
#' @param main.seed An integer to set seed for clustering. This parameter ensures reproduciblity and can be
#'                  kept as default.
#' @param ... Additional arguments to be supplied to clusterNet. *Please see details for information regarding additional
#'             required parameters for each type*
#' @details
#' ## Additional Parameters
#' **weighted** *REQUIRED only for type = "adjacency_matrix"* TRUE/FALSE indicating whether the adjacency graph is weighted.
#'
#' **metaboliteA** *REQUIRED only for type = "edge_list"* A string indicating the column name of the first column in the edge pair.
#'
#' **metaboliteB** *REQUIRED only for type = "edge_list"* A string indicating the column name of the second column in the edge pair.
#'
#'
#' @return A list containing:
#'
#'         1. A data frame, named "subnetworks", containing sub-network determinations for the nodes within the input network.
#'
#'         2. A data frame, named "summary", with a summary of the number of edges and nodes in each subnetwork.
#'
#' @import igraph
#' @include internals.R
#' @include utilities.R
#' @export
clusterNet <- function(data,
                       type = c("adjacency_matrix", "edge_list"),
                       tau = 0.5,
                       main.seed = 417,
                       ...){


  ###############################
  #**  tau must be above 0.5  **#
  ###############################
  if(tau < 0.5 | tau > 1.0) stop(paste('tau must be greater than 0.5!,
                                       Clustering results below this threshold are not reliable -',
                                       'Please see user documentation for more information!'))

  ##############################
  #**. grab optional params  **#
  ##############################
  #type <- match.arg(type)
  edgeListParams <- list(...)
  message("Clustering will be performed in: UN-DIRECTED mode.")

  ##############################
  #**  create igraph object  **#
  ##############################

  if(type == "adjacency_matrix"){

    #warn if any unused params
    if(!all(names(edgeListParams) %in% c("weighted"))){

      warning(paste0("The following unused parameters were discarded: ",
                     paste(names(edgeListParams)[!(names(edgeListParams) %in%
                                                     c("weighted"))],
                           collapse = ", ")))

    }

    #create graph from adjacency matrix input
    inputNetwork <- graph_from_adjacency_matrix(as.matrix(data), mode = "undirected",
                                                weighted = edgeListParams$weighted)
    metabNames <- V(inputNetwork)$name

  } else if(type == "edge_list"){

    #warn if any unused params
    if(!all(names(edgeListParams) %in% c("metaboliteA", "metaboliteB"))){

      warning(paste0("The following unused parameters were discarded: ",
                     paste(names(edgeListParams)[!(names(edgeListParams) %in%
                                                     c("metaboliteA", "metaboliteB"))],
                           collapse = ", ")))
    }

    #create graph from edge list input
    inputNetwork <- graph_from_edgelist(as.matrix(data[,c(edgeListParams$metaboliteA, edgeListParams$metaboliteB)]),
                                        directed = FALSE)
    metabNames <- V(inputNetwork)$name

  }

  #modify the graph
  E(inputNetwork)$lty <- 1
  E(inputNetwork)$color <- "black"

  #more changes if in weighted mode
  if(isTRUE(edgeListParams$weighted) & type == "adjacency_matrix") {

    E(inputNetwork)$lty[is.na(E(inputNetwork)$weight)] <- 2
    E(inputNetwork)$color[is.na(E(inputNetwork)$weight)] <- "red"

  }
  #######################################
  #** run consensus cluster algorithm **#
  #######################################
  fit <- run_consensus_cluster(inputNetwork, tau=tau, main.seed = main.seed)
  consensus_membership <- fit$dcl

  #gather results
  subnetwork_results <- gatherConsensusMatrix(consensus_membership)

  #####################################
  #**Concatenate results for output **#
  #####################################

  summary_list <- list()
  for (loop_cluster in 1:nrow(subnetwork_results) ){
    cluster_c <- induced.subgraph(inputNetwork, V(inputNetwork)$name[(subnetwork_results[loop_cluster,]==1)])
    summary_list[[loop_cluster]] <- data.frame("number_of_nodes"=length(V(cluster_c)),
                                               "number_of_edges"=length(E(cluster_c)),check.names = FALSE)
  }

  #compile summary table
  summary_stat <- data.frame("Subnetworks"= rownames(subnetwork_results), do.call(rbind, summary_list), check.names = FALSE)

  #add names to consensus_membership
  names(consensus_membership) <- metabNames

  #compile output list
  output <- list()
  output$subnetworks <- consensus_membership
  output$summary <- summary_stat
  output <- newClusterNetResults(output)
  return(output)

}

