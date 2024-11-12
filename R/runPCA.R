#' Perform Principal components analysis with `scrapper`
#'
#' @param object A matrix-like object. Rows are features and columns are cells.
#' @export
runPCA <- function(object, ...) UseMethod("runPCA")

#' @param subset_row Integer, logical or character vector specifying which
#' features to use in the PCA (e.g., highly variable genes). If `NULL`, all
#' features in `sce_raw` are used.
#' @param d Integer scalar specifying the number of top PCs to obtain.
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
runPCA.default <- function(object,
                           d = 50L, scale = FALSE, subset_row = NULL, ...,
                           block = NULL, block_weight_policy = NULL,
                           variable_block_weight = c(0, 1000),
                           from_residuals = FALSE, extra_work = 7,
                           iterations = 1000, seed = NULL,
                           realized = TRUE, threads = NULL) {
    rlang::check_dots_empty()
    threads <- set_threads(threads)
    if (!is.null(subset_row)) object <- object[subset_row, , drop = FALSE]
    seed <- check_seed(seed)
    pcs <- scrapper::runPca(
        x = object,
        number = d,
        scale = scale,
        block = block,
        block.weight.policy = block_weight_policy,
        variable.block.weight = variable_block_weight,
        components.from.residuals = from_residuals,
        extra.work = extra_work,
        iterations = iterations,
        seed = seed,
        realized = realized,
        num.threads = threads
    )

    # batch_pcs$components:
    # a numeric matrix containing the top principal components. Each row
    # corresponds to a PC and each column corresponds to a cell.
    out <- t(pcs$components)
    attr(out, "percentVar") <- pcs$variance.explained / pcs$total.variance
    attr(out, "totalVariance") <- pcs$total.variance
    attr(out, "rotation") <- pcs$rotation
    out
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
