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

#' @param tips The cell `index`/`indices` from which to calculate the `DPT(s)`
#' (integer of length `[1, 3]`).
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

#' Find roots in a DiffusionMap object
#'
#' @description
#' Find roots in a DiffusionMap object. If we know the `start` and `end` of a
#' root, we then check if the tips of the root are in the start and ends node we
#' already known (see [`find_tips`][destiny::find_tips]). This function in this
#' way check `n_root` nodes with the top `DPT` in the start cluster and return
#' roots meet this criteria.
#'
#' @param object A [`DiffusionMap`][destiny::DiffusionMap] object.
#' @param start,ends The `start` and `ends` cluster identity. start must have a
#'   length `1L`, while the length of ends ranges from `1L` to `2L`. All `start`
#'   and `ends` must exist in `ref`. If `ends` is `NULL`, this only check if one
#'   tip is in the start cluster.
#' @param ref All cell identity of `object`, this must have the same length of
#'   `object@@d`.
#' @param n_root The number of nodes we should test.
#' @return An integer index.
#' @export
diffusionMapRoots <- function(object, ...) {
    rlang::check_installed("destiny", "to use `diffusionMapRoots()`")
    UseMethod("diffusionMapRoots")
}

#' @export
diffusionMapRoots.DiffusionMap <- function(object, start, ends = NULL,
                                           ref, n_root = 100L) {
    assert_number_whole(start, min = 1)
    assert_number_whole(n_root, min = 1)
    if (length(object@d) != length(ref)) {
        cli::cli_abort(
            "{.arg ref} must have the same length {.code length(object@d)} ({length(object@d)})"
        )
    }
    if (!all(c(start, ends) %in% ref)) {
        cli::cli_abort("{.arg start} and {.arg ends} must exist in {.arg ref}.")
    }
    if (!is.null(ends)) {
        if (length(ends) > 2L || length(ends) < 1L) {
            cli::cli_abort("{.arg ends} must be of length [1, 2]")
        }
    }

    # we extract the top `n_root` with maximal DPT from start notes.
    roots <- which(ref == start, useNames = FALSE)
    dpt <- methods::new("DPT",
        branch = matrix(), tips = matrix(),
        dm = object
    )
    dpt <- dpt[sample(roots, size = 1L)][roots]
    n_root <- min(n_root, length(roots), na.rm = TRUE)
    roots <- roots[order(dpt, decreasing = TRUE)][seq_len(n_root)]

    roots[vapply(
        roots, is_matched_root, logical(1L),
        dm = object, start = start,
        ends = ends, ref = ref
    )]
}

# root is the index
# start, end1 and end2 are all a single string.
# ref is the reference of start, end1 and end2.
is_matched_root <- function(dm, root, start, ends, ref) {
    tips <- destiny::find_tips(dm, root = root)
    if (is.null(ends)) {
        any(tips %in% start, na.rm = TRUE)
    } else {
        setequal(ref[tips], c(start, ends))
    }
}
