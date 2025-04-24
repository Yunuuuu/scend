#' @inherit scrapper::computeRnaQcMetrics title description
#' @seealso [`computeRnaQcMetrics`][scrapper::computeRnaQcMetrics]
#' @return A data frame with the columns `sum`, `detected`, `subsets_*`, each of
#' which is a numeric vector containing the QC metrics for all cells.
#' @export
perCellQCMetrics <- function(object, ...) {
    UseMethod("perCellQCMetrics")
}

#' @inheritParams modelGeneVar
#' @inheritParams scrapper::computeRnaQcMetrics
#' @export
#' @rdname perCellQCMetrics
perCellQCMetrics.default <- function(object, subsets = NULL, ...,
                                     threads = NULL) {
    rlang::check_dots_empty()
    threads <- set_threads(threads)
    metrics <- scrapper::computeRnaQcMetrics(object,
        subsets = subsets %||% list(), num.threads = threads
    )
    out <- .subset(metrics, c("sum", "detected"))
    if (length(subsets)) {
        subsets <- .subset2(metrics, "subsets")
        names(subsets) <- paste("subsets", names(subsets), sep = "_")
        out <- c(out, subsets)
    }
    quickdf(out)
}

#' @export
perCellQCMetrics.HDF5Matrix <- function(object, ..., threads = NULL) {
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
#' @rdname perCellQCMetrics
perCellQCMetrics.SingleCellExperiment <- function(object, ...,
                                                  assay = "counts") {
    mat <- .get_mat_from_sce(object, assay)
    perCellQCMetrics(object = mat, ...)
}

#' @export
#' @rdname perCellQCMetrics
perCellQCMetrics.Seurat <- function(object, ..., assay = NULL, layer = "data") {
    mat <- .get_mat_from_seurat(object, assay, layer)
    perCellQCMetrics(object = mat, ...)
}

#############################################################
#' Add QC metrics
#' @inheritParams modelGeneVar
#' @inheritDotParams perCellQCMetrics
#' @seealso
#' - [`perCellQCMetrics`]
#' - [`computeRnaQcMetrics`][scrapper::computeRnaQcMetrics]
#' @export
addPerCellQCMetrics <- function(object, ...) UseMethod("addPerCellQCMetrics")

#' @export
addPerCellQCMetrics.SingleCellExperiment <- function(object, ...) {
    out <- perCellQCMetrics(object = object, ...)
    SummarizedExperiment::colData(object) <- cbind(
        SummarizedExperiment::colData(object), out
    )
    object
}

#' @export
addPerCellQCMetrics.Seurat <- function(object, ...) {
    out <- perCellQCMetrics(object = object, ...)
    SeuratObject::AddMetaData(object, out)
}
