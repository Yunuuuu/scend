#' Quick clustering of cells with `scrapper`
#'
#' @inheritDotParams clusterSNNGraph
#' @inheritParams runPCA
#' @inheritParams logNormCounts
#' @param n_hvgs Integer specifying the number of top genes to retain.
#' @param prop_hvgs Numeric scalar specifying the proportion of genes to report
#' as HVGs.
#' @param seed Integer specifying the seed to use for the initial random vector
#' in `IRLBA` and for `multi-level` or `Leiden` clustering.
#' @seealso
#' - [`logNormCounts`]
#' - [`runPCA`]
#' - [`clusterSNNGraph`]
#' @export
quickCluster <- function(object, ...) UseMethod("quickCluster")

.quickCluster <- function(object, size_factors = NULL,
                          # `runPCA` arguments
                          n_dim = 50L, scale = FALSE,
                          subset_row = NULL, n_hvgs = 2000, prop_hvgs = NULL,
                          from_residuals = FALSE, extra_work = 7,
                          iterations = 1000,
                          # `clusterSNNGraph` arguments
                          ..., seed = NULL, realized = TRUE, threads = NULL) {
    # nromalization, adjust for differences in sequencing depth
    object <- logNormCounts(
        object = object, threads = threads,
        size_factors = size_factors
    )
    seed <- check_seed(seed, 2L)

    if (is.null(subset_row)) {
        fit <- modelGeneVar(object, threads = threads)
        # At least 500 genes, or 10% of genes; whichever is larger.
        subset_row <- getTopHVGs(fit, n = n_hvgs, prop = prop_hvgs)
    }

    # dimensionality reduction
    object <- runPCA(
        object = object,
        threads = threads,
        n_dim = n_dim, scale = scale, subset_row = subset_row,
        from_residuals = from_residuals, extra_work = extra_work,
        iterations = iterations, seed = seed[1L],
        realized = realized
    )

    # run clusterSNNGraph
    clusterSNNGraph(
        object = object, dimred = "PCA",
        ..., threads = threads, seed = seed[2L]
    )
}

#' @export
#' @rdname quickCluster
quickCluster.SingleCellExperiment <- .quickCluster

#' @export
#' @rdname quickCluster
quickCluster.Seurat <- .quickCluster
