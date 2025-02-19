#' Create a diffusion map of cells
#'
#' @inheritParams runPCA
#' @inheritParams destiny::DPT
#' @inheritDotParams destiny::DiffusionMap -data -n_eigs -rotate -n_pcs -suppress_dpt
#' @return
#'  - `default` method: A numeric matrix where rows are cells and columns are
#'                     the two dimensions of the embedding.
#'  - `SingleCellExperiment` method: embedding was added into
#'    [`reducedDims`][SingleCellExperiment::reducedDim] named as `name`.
#'  - `Seurat` method: embedding was added into
#'    [`reductions`][SeuratObject::Seurat-class] named as `name`.
#'
#' In both cases, the attributes of the `DiffusionMap` coordinate matrix contain
#' the `DPT`: the raw output of [`DPT`][destiny::DPT].
#'
#' @export
runDiffusionMap <- function(object, ...) {
    rlang::check_installed("destiny", "to use `runDiffusionMap()`")
    UseMethod("runDiffusionMap")
}

#' @export
#' @rdname runDiffusionMap
runDiffusionMap.default <- function(object, d = 20L, ...,
                                    tips = NULL, w_width = 0.1) {
    dm <- destiny::DiffusionMap(
        data = object, ...,
        n_eigs = d, rotate = FALSE,
        n_pcs = NA, suppress_dpt = FALSE
    )
    if (is.null(tips)) {
        dpt <- destiny::DPT(dm, w_width = w_width)
    } else {
        dpt <- destiny::DPT(dm, tips = tips, w_width = w_width)
    }
    attr(dpt, "w_width") <- w_width
    evs <- destiny::eigenvectors(dm)
    SingleCellExperiment::reduced.dim.matrix(evs, DPT = dpt)
}

#' @export
#' @rdname runDiffusionMap
runDiffusionMap.SingleCellExperiment <- function(object, ...,
                                                 assay = "logcounts",
                                                 dimred = NULL, n_dimred = NULL,
                                                 name = "DiffusionMap") {
    mat <- .get_mat_from_sce(object, assay, dimred, n_dimred)
    dm <- runDiffusionMap(object = mat, ...)
    add_dimred_to_sce(object, dm, name)
}

#' @export
#' @rdname runDiffusionMap
runDiffusionMap.Seurat <- function(object, ...,
                                   assay = NULL, layer = "data",
                                   dimred = NULL, n_dimred = NULL,
                                   name = "DiffusionMap") {
    mat <- .get_mat_from_seurat(object, assay, layer, dimred, n_dimred)
    dm <- runDiffusionMap(object = mat, ...)
    add_dimred_to_seurat(object, dm, name, assay, layer, dimred)
}
