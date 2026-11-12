#' Compute the t-stochastic neighbor embedding with `scrapper`
#'
#' @export
runTSNE <- function(object, ...) {
    UseMethod("runTSNE")
}

#' @param perplexity Numeric scalar specifying the perplexity to use in the
#' t-SNE algorithm.
#' @param n_neighbors Integer scalar specifying the number of neighbors,
#' passed to scrapper. If `NULL`, scrapper derives it from `perplexity`.
#' @param max_depth Integer scalar specifying the maximum depth of the
#' Barnes-Hut quadtree. Smaller values (7-10) improve speed at the cost of
#' accuracy.
#' @param approximate Logical scalar indicating whether to use the
#' \dQuote{leaf approximation} approach, which sacrifices some accuracy for
#' greater speed. Only effective when `max_depth` is small enough for multiple
#' cells to be assigned to the same leaf node of the quadtree.
#' @param max_iter Integer scalar specifying the maximum number of iterations to
#' perform.
#' @param bnparam A \linkS4class{BiocNeighborParam} object specifying the
#' algorithm to use.
#' @param seed Integer seed. If `NULL`, a seed is drawn from R's random number
#' generator, so results can be reproduced with `set.seed()`.
#' @inheritParams modelGeneVar
#' @inheritParams scrapper::runTsne
#' @inherit runPCA return
#' @seealso [`runTsne`][scrapper::runTsne]
#' @importFrom BiocNeighbors AnnoyParam
#' @export
#' @rdname runTSNE
runTSNE.default <- function(object, perplexity = 30L, n_neighbors = NULL,
                            ...,
                            max_depth = 20L, approximate = FALSE,
                            max_iter = 500L, seed = NULL,
                            bnparam = AnnoyParam(), threads = NULL) {
    rlang::check_dots_empty()
    threads <- set_threads(threads)
    seed <- check_seed(seed)
    args <- list(
        x = object,
        perplexity = perplexity,
        max.depth = max_depth,
        leaf.approximation = approximate,
        max.iterations = max_iter,
        seed = seed, num.threads = threads,
        BNPARAM = bnparam
    )
    if (!is.null(n_neighbors)) {
        args$num.neighbors <- n_neighbors
    }
    rlang::inject(scrapper::runTsne(!!!args))
}

#' @export
runTSNE.HDF5Matrix <- function(object, ..., threads = NULL) {
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
#' @rdname runTSNE
runTSNE.SingleCellExperiment <- function(object, ...,
                                         dimred = "PCA", n_dimred = NULL,
                                         assay = NULL,
                                         name = "TSNE") {
    mat <- .get_mat_from_sce(object, assay, dimred, n_dimred)
    tsne <- runTSNE(object = mat, ...)
    add_dimred_to_sce(object, tsne, name)
}

#' @inheritParams runPCA
#' @export
#' @rdname runTSNE
runTSNE.Seurat <- function(object, ...,
                           dimred = "PCA", n_dimred = NULL,
                           assay = NULL, layer = NULL,
                           name = "TSNE") {
    mat <- .get_mat_from_seurat(object, assay, layer, dimred, n_dimred)
    tsne <- runTSNE(object = mat, ...)
    add_dimred_to_seurat(object, tsne, name, assay, layer, dimred)
}
