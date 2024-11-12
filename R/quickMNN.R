#' Quick fastMNN with `scrapper`
#'
#' @inheritDotParams runPCA
#' @inheritParams runMNN
#' @inheritParams logNormCounts
#' @seealso
#' - [`logNormCounts`]
#' - [`runPCA`]
#' - [`runMNN`]
#' @export
quickMNN <- function(object, ...) UseMethod("quickMNN")

#' @export
#' @rdname quickMNN
quickMNN.SingleCellExperiment <- function(
    object, block, size_factors = NULL,
    # `runPCA` arguments
    ...,
    # `runMNN` arguments
    k = 15L, n_mads = 3L,
    robust_iterations = 2, robust_trim = 0.25,
    mass_cap = NULL,
    order = NULL, reference_policy = NULL,
    BNPARAM = AnnoyParam(), threads = NULL,
    name = "corrected") {
    # nromalization, adjust for differences in sequencing depth
    object <- logNormCounts(
        object = object, block = block, threads = threads,
        size_factors = size_factors,
        name = "multiBatchNorm"
    )

    # dimensionality reduction
    object <- runPCA(
        object = object, block = block, ...,
        assay = "multiBatchNorm",
        threads = threads, name = "PCA"
    )

    # run MNN for batch correction
    runMNN(
        object = object, block = block, dimred = "PCA",
        k = k, n_mads = n_mads,
        robust_iterations = robust_iterations,
        robust_trim = robust_trim,
        mass_cap = mass_cap,
        order = order, reference_policy = reference_policy,
        BNPARAM = BNPARAM, threads = threads, name = name
    )
}

#' @export
#' @rdname quickMNN
quickMNN.Seurat <- function(
    object, block, size_factors = NULL,
    # `runPCA` arguments
    ...,
    # `runMNN` arguments
    k = 15L, n_mads = 3L,
    robust_iterations = 2, robust_trim = 0.25,
    mass_cap = NULL,
    order = NULL, reference_policy = NULL,
    BNPARAM = AnnoyParam(), threads = NULL,
    name = "corrected") {
    # nromalization, adjust for differences in sequencing depth
    object <- logNormCounts(
        object = object, block = block, threads = threads,
        size_factors = size_factors,
        name = "multiBatchNorm"
    )

    # dimensionality reduction
    object <- runPCA(
        object = object, block = block, ...,
        layer = "multiBatchNorm",
        threads = threads, name = "PCA"
    )

    # run MNN for batch correction
    runMNN(
        object = object, block = block, dimred = "PCA",
        k = k, n_mads = n_mads,
        robust_iterations = robust_iterations,
        robust_trim = robust_trim,
        mass_cap = mass_cap,
        order = order, reference_policy = reference_policy,
        BNPARAM = BNPARAM, threads = threads, name = name
    )
}
