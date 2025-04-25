#' Fast mutual nearest neighbors correction with `scrapper`
#'
#' @export
runMNN <- function(object, ...) {
    UseMethod("runMNN")
}

#' @param block Factor specifying the block of origin (e.g., batch, sample) for
#' each cell in `object`.
#' @param n_neighbors Integer scalar specifying the number of neighbors to use
#' when identifying MNN pairs.
#' @param n_mads Numeric scalar specifying the number of median absolute
#' deviations to use for removing outliers in the center-of-mass calculations.
#' @param robust_iterations Integer scalar specifying the number of iterations
#' for robust calculation of the center of mass.
#' @param robust_trim Numeric scalar in [0, 1) specifying the trimming
#' proportion for robust calculation of the center of mass.
#' @param mass_cap Integer scalar specifying the cap on the number of
#' observations to use for center-of-mass calculations on the reference dataset.
#' A value of `100,000` may be appropriate for speeding up correction of very
#' large datasets. If `NULL`, no cap is used.
#' @param order Vector containing levels of `block` in the desired merge
#' order. If `NULL`, a suitable merge order is automatically determined.
#' @param reference_policy String specifying the policy to use to choose the
#' first reference batch. This can be based on the largest batch (`"max-size"`),
#' the most variable batch (`"max-variance"`), the batch with the largest
#' residual sum of squares (`"max-rss"`), or the first specified input
#' (`"input"`). Only used for automatic merges, i.e., when `order=NULL`.
#' @inheritParams scrapper::correctMnn
#' @inheritParams runPCA
#' @inheritParams runTSNE
#' @seealso [`correctMnn`][scrapper::correctMnn]
#' @importFrom BiocNeighbors AnnoyParam
#' @export
#' @rdname runMNN
runMNN.default <- function(object, block, n_neighbors = 15L, ...,
                           n_mads = 3L,
                           robust_iterations = 2, robust_trim = 0.25,
                           mass_cap = NULL,
                           order = NULL, reference_policy = NULL,
                           bnparam = AnnoyParam(), threads = NULL) {
    rlang::check_dots_empty()
    # run MNN --------------------------------------------------------
    .runMNN(
        object = object,
        block = block, n_neighbors = n_neighbors,
        n_mads = n_mads,
        robust_iterations = robust_iterations,
        robust_trim = robust_trim,
        mass_cap = mass_cap,
        order = order,
        reference_policy = reference_policy,
        bnparam = bnparam,
        threads = set_threads(threads)
    )
}

#' @export
runMNN.HDF5Matrix <- function(object, ..., threads = NULL) {
    threads <- set_threads(threads)
    if (threads > 1L) {
        cli::cli_warn(c(
            "Cannot use multiple threads for {.cls HDF5Matrix}",
            i = "Will use {.code threads = 1} instead"
        ))
        threads <- 1L
    }
    NextMethod()
}

#' @export
#' @rdname runMNN
runMNN.SingleCellExperiment <- function(object, ...,
                                        dimred = "PCA", n_dimred = NULL,
                                        assay = NULL,
                                        name = "corrected") {
    mat <- .get_mat_from_sce(object, assay, dimred, n_dimred)
    mnn <- runMNN(object = mat, ...)
    add_dimred_to_sce(object, mnn, name)
}

#' @export
#' @rdname runMNN
runMNN.Seurat <- function(object, ...,
                          dimred = "PCA", n_dimred = NULL,
                          assay = NULL, layer = NULL,
                          name = "corrected") {
    mat <- .get_mat_from_seurat(object, assay, layer, dimred, n_dimred)
    mnn <- runMNN(object = mat, ...)
    add_dimred_to_seurat(object, mnn, name, assay, layer, dimred)
}

.runMNN <- function(object, block, n_neighbors, n_mads,
                    robust_iterations, robust_trim, mass_cap,
                    order, reference_policy, bnparam,
                    threads) {
    mnn_res <- scrapper::correctMnn(
        x = object,
        block = block, num.neighbors = n_neighbors,
        num.mads = n_mads,
        robust.iterations = robust_iterations,
        robust.trim = robust_trim,
        mass.cap = mass_cap,
        order = order,
        reference.policy = reference_policy,
        BNPARAM = bnparam,
        num.threads = threads
    )
    new_dimred(
        t(mnn_res$corrected),
        merge_order = mnn_res$merge.order,
        n_pairs = mnn_res$num.pairs
    )
}
