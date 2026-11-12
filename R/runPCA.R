#' Perform Principal components analysis with `scrapper`
#'
#' @param object A matrix-like object. Rows are features and columns are cells.
#' @export
runPCA <- function(object, ...) {
    UseMethod("runPCA")
}

#' @param subset_row Integer, logical or character vector specifying which
#' features to use in the PCA (e.g., highly variable genes). If `NULL`, all
#' features in `object` are used.
#' @param n_dim Integer scalar specifying the number of top PCs to obtain.
#' @param scale Logical scalar indicating whether to scale rows to unit
#' variance.
#' @inheritParams modelGeneVar
#' @inheritParams scrapper::runPca
#' @param from_residuals Logical scalar indicating whether to compute the PC
#' scores from the residuals in the presence of a blocking factor. By default,
#' the residuals are only used to compute the rotation matrix, and the original
#' expression values of the cells are projected onto this new space. Only used
#' if `block` is not `NULL`.
#' @param extra_work Integer scalar specifying the extra dimensions for the
#' IRLBA workspace.
#' @param subset Features to use for PCA, passed to scrapper after applying
#' `subset_row`. Unlike `subset_row`, this retains a rotation for all input
#' features. See [scrapper::runPca].
#' @param tolerance Convergence tolerance for IRLBA. If `NULL`, use scrapper's
#' default. Setting this requires a scrapper version that supports `tolerance`.
#' @param realized Logical scalar indicating whether to realize `object` into an
#' optimal memory layout for `IRLBA`. This improves computation time at the cost
#' of increased memory usage.
#' @param seed Integer seed. If `NULL`, a seed is drawn from R's random number
#' generator, so results can be reproduced with `set.seed()`.
#' @seealso [`runPca`][scrapper::runPca]
#' @return
#'  - `default` method: A numeric matrix where rows are cells and columns are
#'                    the two dimensions of the embedding.
#'  - `SingleCellExperiment` method: embedding was added into
#'    [`reducedDims`][SingleCellExperiment::reducedDim] named as `name`.
#'  - `Seurat` method: embedding was added into
#'    [`reductions`][SeuratObject::Seurat-class] named as `name`.
#' @export
#' @rdname runPCA
runPCA.default <- function(object, n_dim = 50L, scale = FALSE,
                           subset_row = NULL, ...,
                           block = NULL, block_weight_policy = NULL,
                           variable_block_weight = c(0, 1000),
                           from_residuals = FALSE, extra_work = 7,
                           iterations = 1000, seed = NULL,
                           realized = TRUE, subset = NULL,
                           tolerance = NULL, warn = TRUE,
                           threads = NULL) {
    rlang::check_dots_empty()
    threads <- set_threads(threads)
    if (!is.null(subset_row)) object <- object[subset_row, , drop = FALSE]
    seed <- check_seed(seed)
    args <- list(
        x = object,
        number = n_dim,
        scale = scale,
        subset = subset,
        block = block,
        block.weight.policy = block_weight_policy,
        variable.block.weight = variable_block_weight,
        extra.work = extra_work,
        iterations = iterations,
        warn = warn,
        seed = seed,
        realized = realized,
        num.threads = threads
    )
    if (!is.null(tolerance)) args$tolerance <- tolerance
    residual_arg <- if (utils::packageVersion("scrapper") >= "1.7.1") {
        "center.scores.by.block"
    } else {
        "components.from.residuals"
    }
    args[[residual_arg]] <- from_residuals
    pcs <- rlang::inject(scrapper::runPca(!!!args))

    # batch_pcs$components:
    # a numeric matrix containing the top principal components. Each row
    # corresponds to a PC and each column corresponds to a cell.
    new_dimred(
        t(pcs$components),
        percentVar = pcs$variance.explained / pcs$total.variance,
        totalVariance = pcs$total.variance,
        rotation = pcs$rotation
    )
}

#' @export
runPCA.HDF5Matrix <- function(object, ..., threads = NULL) {
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

#' @param dimred String or integer scalar specifying the existing dimensionality
#' reduction results to use.
#' @param n_dimred Integer scalar or vector specifying the dimensions to use if
#' dimred is specified.
#' @param name String specifying the name to be used to store the result in the
#' [`reducedDims`][SingleCellExperiment::reducedDim] or
#' [`reductions`][SeuratObject::Seurat-class] of the output.
#' @export
#' @rdname runPCA
runPCA.SingleCellExperiment <- function(object, ...,
                                        assay = "logcounts",
                                        dimred = NULL, n_dimred = NULL,
                                        name = "PCA") {
    mat <- .get_mat_from_sce(object, assay, dimred, n_dimred)
    pca <- runPCA(object = mat, ...)
    add_dimred_to_sce(object, pca, name)
}

#' @export
#' @rdname runPCA
runPCA.Seurat <- function(object, ...,
                          assay = NULL, layer = "data",
                          dimred = NULL, n_dimred = NULL,
                          name = "PCA") {
    mat <- .get_mat_from_seurat(object, assay, layer, dimred, n_dimred)
    pca <- runPCA(object = mat, ...)
    add_dimred_to_seurat(object, pca, name, assay, layer, dimred)
}
