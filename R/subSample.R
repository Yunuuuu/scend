#' @inherit scrapper::subsampleByNeighbors title description return
#' @inheritParams modelGeneVar
#' @export
subSample <- function(object, ...) {
    UseMethod("subSample")
}

#' @export
#' @rdname subSample
subSample.SingleCellExperiment <- function(object, ..., assay = "counts") {
    mat <- .get_mat_from_sce(object, assay)
    subSample(object = mat, ...)
}

#' @export
#' @rdname subSample
subSample.Seurat <- function(object, ...,
                             assay = NULL, layer = "counts") {
    mat <- .get_mat_from_seurat(object, assay, layer)
    subSample(object = mat, ...)
}

#' @param n_neighbors Integer scalar specifying the number of neighbors to
#' use. Larger values result in greater downsampling. Only used if `object` does
#' not contain existing nearest-neighbor results.
#' @param min_remaining Integer scalar specifying the minimum number of
#' remaining (i.e., unselected) neighbors that a cell must have in order to be
#' considered for selection. This should be less than or equal to
#' `n_neighbors`.
#' @inheritParams runTSNE
#' @seealso [subsampleByNeighbors][scrapper::subsampleByNeighbors]
#' @importFrom BiocNeighbors AnnoyParam
#' @export
#' @rdname subSample
subSample.default <- function(object, n_neighbors = 20, ...,
                              min_remaining = 10, BNPARAM = AnnoyParam(),
                              threads = NULL) {
    rlang::check_dots_empty()
    assert_number_whole(n_neighbors)
    assert_number_whole(min_remaining)
    threads <- set_threads(threads)
    scrapper::subsampleByNeighbors(
        x = object,
        num.neighbors = n_neighbors,
        min.remaining = min_remaining,
        num.threads = threads,
        BNPARAM = BNPARAM
    )
}

#' @export
subSample.HDF5Matrix <- function(object, ..., threads = NULL) {
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
