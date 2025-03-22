#' Compute library size factors
#'
#' @inheritParams modelGeneVar
#' @export
librarySizeFactors <- function(object, ...) {
    UseMethod("librarySizeFactors")
}

#' @export
#' @rdname librarySizeFactors
librarySizeFactors.default <- function(object, ..., threads = NULL) {
    rlang::check_dots_empty()
    threads <- set_threads(threads)
    # `scuttle::librarySizeFactors()`
    size_factors <- scrapper::computeRnaQcMetrics(object,
        subsets = list(), num.threads = threads
    )
    .subset2(size_factors, "sum")
}

#' @export
#' @rdname librarySizeFactors
librarySizeFactors.DelayedArray <- function(object, ...) {
    rlang::check_dots_empty()
    # Will fallback to `DelayedArray::blockApply()`
    colSums(object)
}

#' @export
#' @rdname librarySizeFactors
librarySizeFactors.SummarizedExperiment <- function(object,
                                                    ...,
                                                    assay = "counts") {
    librarySizeFactors(object = .get_mat_from_sce(object, assay), ...)
}

#' @export
#' @rdname librarySizeFactors
librarySizeFactors.Seurat <- function(object, ..., assay = NULL,
                                      layer = "counts") {
    mat <- .get_mat_from_seurat(object, assay, layer)
    librarySizeFactors(object = mat, ...)
}
