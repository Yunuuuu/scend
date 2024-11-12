#' @inherit scrapper::normalizeCounts title description
#' @export
logNormCounts <- function(object, ...) {
    rlang::check_dots_used()
    UseMethod("logNormCounts")
}

#' @inheritParams modelGeneVar
#' @param size_factors A numeric vector of length equal to the number of cells
#' in `object`, containing positive size factors for all cells.
#' @param log Numeric scalar specifying the base of the log-transformation. If
#' `NULL`, won't do log-transformation.
#' @param pseudo_count Numeric scalar specifying the positive pseudo-count to
#' add before any log-transformation. Ignored if `log = NULL`.
#' @param center_size_factors Logical scalar indicating whether size factors
#' should be centered at unity before being used.
#' @param block Vector or factor of length equal to the number of cells,
#' specifying the batch of origin for each cell. Alternatively `NULL` if all
#' cells belong to the same batch. Ignored if `center_size_factors = FALSE`.
#' @param mode String specifying how to scale size factors across blocks.
#' `"lowest"` will scale all size factors by the the lowest per-block average.
#' `"per-block"` will center the size factors in each block separately. Ignored
#' if `center_size_factors = FALSE`.
#' @seealso [normalizeCounts][scrapper::normalizeCounts]
#' @export
#' @rdname logNormCounts
logNormCounts.default <- function(object, size_factors = NULL,
                                  ...,
                                  log = 2L, pseudo_count = 1L,
                                  center_size_factors = TRUE,
                                  block = NULL, mode = NULL,
                                  threads = NULL) {
    ans <- .logNormCounts(
        object = object,
        size_factors = size_factors,
        log = log, pseudo_count = pseudo_count,
        center_size_factors = center_size_factors,
        block = block, mode = mode, threads = threads
    )
    .subset2(ans, "logcounts")
}

#' @importFrom rlang caller_call
.logNormCounts <- function(object, size_factors = NULL,
                           log = 2L, pseudo_count = 1L,
                           center_size_factors = TRUE,
                           block = NULL, mode = NULL,
                           threads = NULL) {
    call <- caller_call()

    assert_(size_factors, is.numeric, "a numeric",
        allow_null = TRUE, call = call
    )
    assert_number_decimal(log, allow_null = TRUE, call = call)
    assert_number_whole(pseudo_count, call = call)
    assert_bool(center_size_factors, call = call)
    if (is.null(mode)) {
        mode <- "per-block"
    } else {
        mode <- rlang::arg_match0(
            mode,
            c("per-block", "lowest"),
            error_call = call
        )
    }
    threads <- set_threads(threads, call = call)
    if (is.null(size_factors)) {
        # `scuttle::librarySizeFactors()`
        size_factors <- scrapper::computeRnaQcMetrics(object,
            subsets = list(), num.threads = threads
        )
        size_factors <- .subset2(size_factors, "sum")
    }
    if (center_size_factors) {
        size_factors <- scrapper::centerSizeFactors(
            size.factors = size_factors, block = block, mode = mode
        )
    }
    logcounts <- scrapper::normalizeCounts(
        x = object,
        size.factors = size_factors,
        log = !is.null(log), pseudo.count = pseudo_count,
        log.base = log, preserve.sparsity = FALSE
    )
    list(logcounts = logcounts, size_factors = size_factors)
}

#' @param assay Integer scalar or string indicating which assay of x
#' contains the expression values.
#' @param name String specifying the name to be used to store the assay name in
#' the [`assays`][SummarizedExperiment::assays] or the layer name in
#' [`SetAssayData`][SeuratObject::SetAssayData].
#' @export
#' @rdname logNormCounts
logNormCounts.SingleCellExperiment <- function(object, size_factors = NULL,
                                               ...,
                                               assay = "counts",
                                               name = "logcounts") {
    internal_sf <- SingleCellExperiment::sizeFactors(object)
    assign_sf <- is.null(size_factors) && is.null(internal_sf)
    size_factors <- size_factors %||% internal_sf
    mat <- .get_mat_from_sce(object, assay)
    ans <- .logNormCounts(object = mat, size_factors = size_factors, ...)
    SummarizedExperiment::assay(object, name) <- .subset2(ans, "logcounts")
    if (assign_sf) {
        SingleCellExperiment::sizeFactors(object) <- .subset2(
            ans, "size_factors"
        )
    }
    object
}

#' @param layer Name of the layer to get from the assay data.
#' @export
#' @rdname logNormCounts
logNormCounts.Seurat <- function(object, ..., assay = NULL,
                                 layer = "counts", name = "data") {
    mat <- .get_mat_from_seurat(object, assay, layer)
    ans <- .logNormCounts(object = mat, ...)
    SeuratObject::SetAssayData(object,
        layer = name, assay = assay,
        new.data = .subset2(ans, "logcounts")
    )
}
