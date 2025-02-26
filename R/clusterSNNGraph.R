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

#' @param ... Additional arguments passed on to
#' [`cluster_louvain`][igraph::cluster_louvain],
#' [`cluster_leiden`][igraph::cluster_leiden] or
#' [`cluster_walktrap`][igraph::cluster_walktrap].
#' @param n_neighbors Integer scalar specifying the number of neighbors to use
#' to construct the graph.
#' @param method String specifying the community detection method to use.
#' Options are multi-level (`"multilevel"`/`"louvain"`), Walktrap (`"walktrap"`)
#' or Leiden (`"leiden"`) (default).
#' @param scheme String specifying the weighting scheme to use for constructing
#' the SNN graph. This can be `"ranked"` (default), `"jaccard"` or `"number"`.
#' @param resolution Numeric scalar specifying the resolution to use for
#' `multi-level` or `Leiden` clustering.
#' @param steps Integer scalar specifying the number of steps to use for
#' `Walktrap` clustering.
#' @param objective String specifying the objective function to use for Leiden
#' clustering: "cpm" or "modularity" (default).
#' @param seed Integer scalar specifying the seed to use for `multi-level` or
#' `Leiden` clustering.
#' @inheritParams scrapper::clusterGraph
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
                                    ...,
                                    scheme = NULL, resolution = 1L,
                                    objective = NULL, steps = 4L,
                                    seed = NULL,
                                    BNPARAM = AnnoyParam(), threads = NULL) {
    threads <- set_threads(threads)
    assert_number_whole(n_neighbors)
    method <- arg_match(
        method,
        c("multilevel", "louvain", "walktrap", "leiden"),
        default = "leiden"
    )
    scheme <- arg_match(scheme, c("ranked", "jaccard", "number"))
    if (identical(method, "walktrap")) {
        assert_number_whole(steps)
    } else {
        assert_number_decimal(resolution)
        if (identical(method, "leiden")) {
            objective <- arg_match(objective, c("modularity", "cpm"))
        }
    }
    seed <- check_seed(seed)
    set_seed(seed)
    graph <- scrapper::buildSnnGraph(
        object,
        num.neighbors = n_neighbors,
        weight.scheme = scheme,
        num.threads = threads,
        BNPARAM = BNPARAM
    )
    weights <- graph$weights
    if (!inherits(graph, "igraph")) {
        my_graph <- igraph::make_undirected_graph(
            graph$edges,
            n = graph$vertices
        )
        igraph::E(my_graph)$weight <- weights
        graph <- my_graph
    }
    if (any(method == c("multilevel", "louvain"))) {
        clustering <- igraph::cluster_louvain(
            graph = graph, resolution = resolution, weights = weights,
            ...
        )
    } else if (method == "leiden") {
        leiden.objective <- if (objective == "cpm") {
            "CPM"
        } else {
            "modularity"
        }
        clustering <- igraph::cluster_leiden(
            graph = graph,
            objective_function = leiden.objective,
            resolution_parameter = resolution,
            weights = weights,
            ...
        )
    } else if (method == "walktrap") {
        clustering <- igraph::cluster_walktrap(
            graph = graph, steps = steps, weights = weights,
            ...
        )
    }
    out <- as.factor(clustering$membership)
    for (i in setdiff(names(clustering), "membership")) {
        attr(out, i) <- .subset2(clustering, i)
    }
    out
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
