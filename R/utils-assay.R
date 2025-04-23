# we always regard features in row and cells in column
.get_mat_from_sce <- function(x, assay, dimred = NULL, n_dimred = NULL) {
    if (!is.null(dimred)) {
        # value is expected to be a matrix or matrix-like object with number of
        # rows equal to ncol(x).
        mat <- SingleCellExperiment::reducedDim(x, dimred)
        if (!is.null(n_dimred)) {
            if (length(n_dimred) == 1L) n_dimred <- seq_len(n_dimred)
            mat <- mat[, n_dimred, drop = FALSE]
        }
        t(mat)
    } else {
        SummarizedExperiment::assay(x, assay)
    }
}

.get_mat_from_seurat <- function(x, assay, layer,
                                 dimred = NULL, n_dimred = NULL) {
    if (!is.null(dimred)) {
        # value is expected to be a matrix or matrix-like object with number of
        # rows equal to ncol(x).
        mat <- SeuratObject::Embeddings(x, reduction = dimred)
        if (!is.null(n_dimred)) {
            if (length(n_dimred) == 1L) n_dimred <- seq_len(n_dimred)
            mat <- mat[, n_dimred, drop = FALSE]
        }
        t(mat)
    } else {
        SeuratObject::GetAssayData(x, assay = assay, layer = layer)
    }
}

.get_assay_from_seurat <- function(x, assay, layer, dimred = NULL) {
    if (!is.null(dimred)) {
        SeuratObject::DefaultAssay(x[[dimred]])
    } else if (!is.null(assay)) {
        assay
    } else {
        SeuratObject::DefaultAssay(x)
    }
}

add_dimred_to_sce <- function(object, value, name) {
    SingleCellExperiment::reducedDim(object, name) <- value
    object
}

add_dimred_to_seurat <- function(object, value, name, assay, layer,
                                 dimred = NULL) {
    reduction_key <- SeuratObject::Key(name, quiet = TRUE)
    object[[name]] <- SeuratObject::CreateDimReducObject(
        embeddings = value,
        stdev = apply(value, 2L, stats::sd, simplify = TRUE),
        assay = .get_assay_from_seurat(object, assay, layer, dimred),
        key = reduction_key
    )
    object
}

new_dimred <- function(dimred, ..., class = NULL) {
    # SingleCellExperiment can only accept `reduced.dim.matrix` class
    class(dimred) <- c(class, "reduced.dim.matrix", "matrix")
    mostattributes(dimred) <- c(attributes(dimred), list(...))
    dimred
}
