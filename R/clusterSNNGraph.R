#' Graph-based clustering with `scrapper`
#'
#' @seealso
#' - [`buildSnnGraph`][scrapper::buildSnnGraph]
#' - [`cluster_louvain`][igraph::cluster_louvain]
#' - [`cluster_leiden`][igraph::cluster_leiden]
#' - [`cluster_walktrap`][igraph::cluster_walktrap]
#' @export
clusterSNNGraph <- function(object, ...) {
    UseMethod("clusterSNNGraph")
}

#' @param ... Additional arguments passed on to [`clusterGraph()`].
#' @inheritParams clusterGraph
#' @param n_neighbors Integer scalar specifying the number of neighbors to use
#' to construct the graph.
#' @param scheme String specifying the weighting scheme to use for constructing
#' the SNN graph. This can be `"ranked"` (default), `"jaccard"` or `"number"`.
#' @inheritParams scrapper::buildSnnGraph
#' @return An integer vector with cluster assignments for each cell. Each method
#' may also return additional attributes.
#' - For method=`"multilevel"`/`"louvain"`, we have:
#'    * `levels`, a list of integer vectors with cluster assignments for each
#' cell at each level. Assignments are sorted by decreasing resolution (i.e.,
#' fewer, larger clusters).
#'    * `modularity`, a numeric vector containing the modularity of each level.
#'    * `best`, the level with the lowest modularity.
#' - For method=`"leiden"`, we have:
#'    * `quality`: a numeric scalar containing the quality of the clustering
#'      (either the modularity or a related score).
#' - For method=`"walktrap"`, we have:
#'    * `merges`: an integer matrix specifying how the clusters were merged to
#'      obtain membership. Each row corresponds to a merge step and contains the
#'      IDs of the temporary clusters (not the same as those in membership).
#'    * `modularity`: a numeric vector containing the modularity before and
#'      after each merge step.
#' @importFrom BiocNeighbors AnnoyParam
#' @export
#' @rdname clusterSNNGraph
clusterSNNGraph.default <- function(object, n_neighbors = 10L, method = NULL,
                                    ..., scheme = NULL,
                                    BNPARAM = AnnoyParam(), threads = NULL) {
    threads <- set_threads(threads)
    assert_number_whole(n_neighbors)
    scheme <- arg_match(scheme, c("ranked", "jaccard", "number"))
    graph <- scrapper::buildSnnGraph(
        object,
        num.neighbors = n_neighbors,
        weight.scheme = scheme,
        num.threads = threads,
        BNPARAM = BNPARAM
    )
    if (!inherits(graph, "igraph")) {
        my_graph <- igraph::make_undirected_graph(
            graph$edges,
            n = graph$vertices
        )
        igraph::E(my_graph)$weight <- graph$weights
        graph <- my_graph
    }
    clusterGraph(object = graph, ..., method = method)
}

#' @inheritParams runPCA
#' @export
#' @rdname clusterSNNGraph
clusterSNNGraph.SingleCellExperiment <- function(object, ...,
                                                 dimred = "PCA",
                                                 n_dimred = NULL,
                                                 assay = NULL) {
    mat <- .get_mat_from_sce(object, assay, dimred, n_dimred)
    clusterSNNGraph(object = mat, ...)
}

#' @inheritParams runPCA
#' @export
#' @rdname clusterSNNGraph
clusterSNNGraph.Seurat <- function(object, ...,
                                   dimred = "PCA", n_dimred = NULL,
                                   assay = NULL, layer = NULL) {
    mat <- .get_mat_from_seurat(object, assay, layer, dimred, n_dimred)
    clusterSNNGraph(object = mat, ...)
}
