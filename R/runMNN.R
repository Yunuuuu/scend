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
#' @param n_steps Integer specifying the number of steps in the recursive
#' neighbor search used to compute the center of mass.
#' @param merge_policy String specifying the merge policy: `"rss"` (default),
#' `"size"`, `"variance"`, or `"input"`. With `"input"`, the last batch is
#' corrected first; arrange the input batches to control the merge order.
#' @inheritParams scrapper::correctMnn
#' @inheritParams runPCA
#' @inheritParams runTSNE
#' @seealso [`correctMnn`][scrapper::correctMnn]
#' @importFrom BiocNeighbors AnnoyParam
#' @export
#' @rdname runMNN
runMNN.default <- function(object, block, n_neighbors = 15L, ...,
                           n_steps = 1L, merge_policy = NULL,
                           bnparam = AnnoyParam(), threads = NULL) {
    rlang::check_dots_empty()
    merge_policy <- arg_match(
        merge_policy,
        c("rss", "size", "variance", "input")
    )
    assert_number_whole(n_steps, min = 1)
    # run MNN --------------------------------------------------------
    .runMNN(
        object = object,
        block = block, n_neighbors = n_neighbors,
        n_steps = n_steps, merge_policy = merge_policy,
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

.runMNN <- function(object, block, n_neighbors, n_steps,
                    merge_policy, bnparam,
                    threads) {
    mnn_res <- scrapper::correctMnn(
        x = object,
        block = block, num.neighbors = n_neighbors,
        num.steps = n_steps,
        merge.policy = merge_policy,
        BNPARAM = bnparam,
        num.threads = threads
    )
    corrected <- t(mnn_res$corrected)
    rownames(corrected) <- colnames(object)
    new_dimred(corrected)
}
