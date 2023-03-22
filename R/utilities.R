#' clusterNetResults
#'
#' clusterNetResults creates a generic for clusterNetResults
#'
#' @title clusterNetResults
#' @param x clusterNetResults object
#' @export
clusterNetResults <- function(x){
  UseMethod("clusterNetResults",x)
}
#' @noRd
#' @export
newClusterNetResults <- function(x = list()){
  stopifnot(is.list(x))
  structure(x, class = "clusterNetResults")
}
#' as.clusterNetResults
#'
#' as.clusterNetResults checks that an object is of class "clusterNetResults"
#' @param x an R variable
#'
#' @return a logical indicating whether the variable is of class "clusterNetResults"
#' @export
as.clusterNetResults <- function(x){
  class(x) <- c("clusterNetResults", class(x))
  x
}
#' custom method to index clusterNetResults
#' @export
`[.clusterNetResults` <- function(x,i){
  newClusterNetResults(NextMethod())
}
