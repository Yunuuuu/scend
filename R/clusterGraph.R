#' Graph-based clustering
#'
#' @param object A [`graph`][igraph::graph] object.
#' @param ... Additional arguments passed on to
#' [`cluster_louvain`][igraph::cluster_louvain],
#' [`cluster_leiden`][igraph::cluster_leiden] or
#' [`cluster_walktrap`][igraph::cluster_walktrap].
#' @param method String specifying the community detection method to use.
#' Options are multi-level (`"multilevel"`/`"louvain"`), Walktrap (`"walktrap"`)
#' or Leiden (`"leiden"`) (default).
#' @param resolution Numeric scalar specifying the resolution to use for
#' `multi-level` or `Leiden` clustering.
#' @param steps Integer scalar specifying the number of steps to use for
#' `Walktrap` clustering.
#' @param objective String specifying the objective function to use for Leiden
#' clustering: "cpm" or "modularity" (default).
#' @param seed Integer scalar specifying the seed to use for `multi-level` or
#' `Leiden` clustering.
#' @return A factor of the cluster assignments. Each method may also return
#' additional attributes.
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
#' @export
clusterGraph <- function(object, ...) UseMethod("clusterGraph")

#' @export
#' @rdname clusterGraph
clusterGraph.igraph <- function(object, method = NULL, ..., resolution = 1L,
                                objective = NULL, steps = 4L, seed = NULL) {
    method <- arg_match(
        method,
        c("multilevel", "louvain", "walktrap", "leiden"),
        default = "leiden"
    )
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
    if (any(method == c("multilevel", "louvain"))) {
        clustering <- igraph::cluster_louvain(
            graph = object, weights = NULL, resolution = resolution,
            ...
        )
    } else if (method == "leiden") {
        leiden.objective <- if (objective == "cpm") {
            "CPM"
        } else {
            "modularity"
        }
        clustering <- igraph::cluster_leiden(
            graph = object,
            objective_function = leiden.objective,
            resolution_parameter = resolution,
            weights = NULL,
            ...
        )
    } else if (method == "walktrap") {
        clustering <- igraph::cluster_walktrap(
            graph = object, weights = NULL, steps = steps,
            ...
        )
    }
    out <- as.factor(clustering$membership)
    for (i in setdiff(names(clustering), "membership")) {
        attr(out, i) <- .subset2(clustering, i)
    }
    out
}
