#' Quick clustering of cells with `scrapper`
#'
#' @inheritDotParams clusterSNNGraph
#' @inheritParams runPCA
#' @inheritParams logNormCounts
#' @param n_hvgs Integer specifying the number of top genes to retain.
#' @param prop_hvgs Numeric scalar specifying the proportion of genes to report
#' as HVGs.
#' @param seed Integer scalar specifying the seed to use for the initial random
#' vector in `IRLBA` and for `multi-level` or `Leiden` clustering.
#' @seealso
#' - [`logNormCounts`]
#' - [`runPCA`]
#' - [`clusterSNNGraph`]
#' @export
quickCluster <- function(object, ...) UseMethod("quickCluster")

.quickCluster <- function(
    object, size_factors = NULL,
    # `runPCA` arguments
    d = 50L, scale = FALSE,
    subset_row = NULL, n_hvgs = 500, prop_hvgs = 0.1,
    block_weight_policy = NULL,
    variable_block_weight = c(0, 1000),
    from_residuals = FALSE, extra_work = 7,
    iterations = 1000, seed = NULL,
    realized = TRUE,
    # `clusterSNNGraph` arguments
    ..., threads = NULL) {
    # nromalization, adjust for differences in sequencing depth
    object <- logNormCounts(
        object = object, threads = threads,
        size_factors = size_factors
    )

    seed <- check_seed(seed)
    set_seed(seed)

    if (is.null(subset_row)) {
        fit <- modelGeneVar(object)
        # At least 500 genes, or 10% of genes; whichever is larger.
        subset_row <- getTopHVGs(fit, n = n_hvgs, prop = prop_hvgs)
    }

    # dimensionality reduction
    object <- runPCA(
        object = object,
        threads = threads,
        d = d, scale = scale, subset_row = subset_row,
        block_weight_policy = block_weight_policy,
        variable_block_weight = variable_block_weight,
        from_residuals = from_residuals, extra_work = extra_work,
        iterations = iterations, seed = random_seed(1L),
        realized = realized
    )

    # run MNN for batch correction
    clusterSNNGraph(
        object = object, dimred = "PCA",
        ..., threads = threads, seed = random_seed(1L)
    )
}

#' @export
#' @rdname quickCluster
quickCluster.SingleCellExperiment <- .quickCluster

#' @export
#' @rdname quickCluster
quickCluster.Seurat <- .quickCluster
