#' Compute the uniform manifold approximation and projection with `scrapper`
#'
#' @export
runUMAP <- function(object, ...) {
    UseMethod("runUMAP")
}

#' @param n_dim Integer scalar specifying the number of dimensions of the output
#' embedding.
#' @param n_neighbors Integer scalar specifying the number of neighbors to use
#' in the UMAP algorithm.
#' @param n_epochs Integer scalar specifying the number of epochs to perform.
#' If `NULL` (or `-1` for compatibility), scrapper chooses an appropriate value.
#' @param min_dist Numeric scalar specifying the minimum distance between
#' points.
#' @param optimization Logical scalar specifying whether to parallelize the
#' optimization step using `threads`. If `FALSE`, optimization uses one thread.
#' @param seed Integer seed, or two integers for initialization and optimization
#' respectively. A single seed is expanded reproducibly into two seeds. If
#' `NULL`, seeds are drawn from R's random number generator.
#' @inheritParams runTSNE
#' @inherit runPCA return
#' @seealso [`runUmap`][scrapper::runUmap]
#' @importFrom BiocNeighbors AnnoyParam
#' @export
#' @rdname runUMAP
runUMAP.default <- function(object, n_dim = 2L, n_neighbors = 15L,
                            n_epochs = NULL, min_dist = 0.01,
                            ..., seed = NULL,
                            bnparam = AnnoyParam(), threads = NULL,
                            optimization = FALSE) {
    rlang::check_dots_empty()
    assert_number_whole(n_epochs, min = -1, allow_null = TRUE)
    assert_bool(optimization)
    if (isTRUE(n_epochs == -1)) n_epochs <- NULL
    threads <- set_threads(threads)
    seed <- check_seed(seed, 2L)
    args <- list(
        x = object,
        num.dim = n_dim,
        num.neighbors = n_neighbors,
        num.epochs = n_epochs,
        min.dist = min_dist,
        initialize.seed = seed[1L],
        optimize.seed = seed[2L],
        num.threads = threads,
        BNPARAM = bnparam
    )
    if (utils::packageVersion("scrapper") >= "1.7.1") {
        args$num.threads.optimize <- if (optimization) threads else 1L
    } else {
        args$parallel.optimization <- optimization
    }
    rlang::inject(scrapper::runUmap(!!!args))
}

#' @export
runUMAP.HDF5Matrix <- function(object, ..., threads = NULL) {
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

#' @inheritParams runPCA
#' @export
#' @rdname runUMAP
runUMAP.SingleCellExperiment <- function(object, ...,
                                         dimred = "PCA", n_dimred = NULL,
                                         assay = NULL,
                                         name = "UMAP") {
    mat <- .get_mat_from_sce(object, assay, dimred, n_dimred)
    umap <- runUMAP(object = mat, ...)
    add_dimred_to_sce(object, umap, name)
}

#' @inheritParams runPCA
#' @export
#' @rdname runUMAP
runUMAP.Seurat <- function(object, ...,
                           dimred = "PCA", n_dimred = NULL,
                           assay = NULL, layer = NULL,
                           name = "UMAP") {
    mat <- .get_mat_from_seurat(object, assay, layer, dimred, n_dimred)
    umap <- runUMAP(object = mat, ...)
    add_dimred_to_seurat(object, umap, name, assay, layer, dimred)
}
