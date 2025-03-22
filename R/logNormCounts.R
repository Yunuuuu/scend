#' @inherit scrapper::normalizeCounts title description
#' @export
logNormCounts <- function(object, ...) {
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
    rlang::check_dots_empty()
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
    threads <- set_threads(threads, call = call)
    if (is.null(size_factors)) {
        # `scuttle::librarySizeFactors()`
        size_factors <- librarySizeFactors(object)
    }
    out <- list(size_factors = size_factors)
    if (center_size_factors) {
        mode <- arg_match(mode, c("per-block", "lowest"), call = call)
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
    c(list(logcounts = logcounts), out)
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
    size_factors <- size_factors %||% SingleCellExperiment::sizeFactors(object)
    mat <- .get_mat_from_sce(object, assay)
    ans <- .logNormCounts(object = mat, size_factors = size_factors, ...)
    SummarizedExperiment::assay(object, name) <- .subset2(ans, "logcounts")
    if (is.null(size_factors)) { # if no size factors
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
