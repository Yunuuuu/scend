#' Compute the uniform manifold approximation and projection with `scrapper`
#'
#' @export
runUMAP <- function(object, ...) UseMethod("runUMAP")

#' @param n_dim Integer scalar specifying the number of dimensions of the output
#' embedding.
#' @param n_neighbors Integer scalar specifying the number of neighbors to use
#' in the UMAP algorithm.
#' @param n_epochs Integer scalar specifying the number of epochs to perform.
#' If set to `-1`, an appropriate number of epochs is chosen based on
#' `ncol(object)`.
#' @param min_dist Numeric scalar specifying the minimum distance between
#' points.
#' @param optimization Logical scalar specifying whether to parallelize the
#' optimization step.
#' @inheritParams runTSNE
#' @inherit runPCA return
#' @seealso [`runUmap`][scrapper::runUmap]
#' @importFrom BiocNeighbors AnnoyParam
#' @export
#' @rdname runUMAP
runUMAP.default <- function(object, n_dim = 2L, n_neighbors = 15L,
                            n_epochs = -1L, min_dist = 0.01,
                            ...,
                            optimization = FALSE, seed = 1234L,
                            BNPARAM = AnnoyParam(), threads = NULL) {
    rlang::check_dots_empty()
    threads <- set_threads(threads)
    scrapper::runUmap(
        x = object,
        num.dim = n_dim,
        num.neighbors = n_neighbors,
        num.epochs = n_epochs,
        min.dist = min_dist,
        seed = seed,
        num.threads = threads,
        parallel.optimization = optimization,
        BNPARAM = BNPARAM
    )
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
