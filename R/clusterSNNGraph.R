#' Graph-based clustering with `scrapper`
#'
#' @seealso
#' - [buildSnnGraph][scrapper::buildSnnGraph]
#' - [clusterGraph][scrapper::clusterGraph]
#' @export
clusterSNNGraph <- function(object, ...) UseMethod("clusterSNNGraph")

#' @param n_neighbors Integer scalar specifying the number of neighbors to use
#' to construct the graph.
#' @param method String specifying the community detection method to use.
#' Options are multi-level (`"multilevel"`), Walktrap (`"walktrap"`) or Leiden
#' (`"leiden"`).
#' @param scheme String specifying the weighting scheme to use for constructing
#' the SNN graph. This can be `"ranked"` (default), `"jaccard"` or `"number"`.
#' @param resolution Numeric scalar specifying the resolution to use for
#' `multi-level` or `Leiden` clustering.
#' @param steps Integer scalar specifying the number of steps to use for
#' `Walktrap` clustering.
#' @param objective String specifying the objective function to use for Leiden
#' clustering: "CPM" or "modularity" (default).
#' @param seed Integer scalar specifying the seed to use for multi-level or
#' Leiden clustering.
#' @inheritParams scrapper::clusterGraph
#' @inheritParams scrapper::buildSnnGraph
#' @return An integer vector with cluster assignments for each cell. Each method
#' may also return additional attributes.
#' - For method=`"multilevel"`, we have:
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
clusterSNNGraph.default <- function(object, n_neighbors = 15L,
                                    method = "leiden",
                                    ...,
                                    scheme = NULL, resolution = 1L,
                                    objective = NULL, steps = 4L,
                                    seed = 1234L,
                                    BNPARAM = AnnoyParam(), threads = NULL) {
    rlang::check_dots_empty()
    threads <- set_threads(threads)
    assert_number_whole(n_neighbors)
    method <- match.arg(method, c("multilevel", "walktrap", "leiden"))
    scheme <- match.arg(scheme, c("ranked", "jaccard", "number"))
    assert_number_decimal(resolution)
    objective <- match.arg(objective, c("modularity", "CPM"))
    assert_number_whole(steps)
    assert_number_whole(seed)
    graph <- scrapper::buildSnnGraph(
        object,
        num.neighbors = n_neighbors,
        weight.scheme = scheme,
        num.threads = threads,
        BNPARAM = BNPARAM
    )
    clustering <- scrapper::clusterGraph(
        x = graph, method = method,
        multilevel.resolution = resolution,
        leiden.resolution = resolution,
        leiden.objective = objective,
        walktrap.steps = steps,
        seed = seed
    )
    out <- clustering$membership
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
