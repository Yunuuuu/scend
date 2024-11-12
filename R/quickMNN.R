#' Quick fastMNN with `scrapper`
#'
#' @inheritDotParams runMNN
#' @inheritParams runPCA
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
    d = 50L, scale = FALSE, subset_row = NULL,
    block_weight_policy = NULL,
    variable_block_weight = c(0, 1000),
    from_residuals = FALSE, extra_work = 7,
    iterations = 1000, seed = NULL,
    realized = TRUE,
    # `runMNN` arguments
    ..., threads = NULL) {
    # nromalization, adjust for differences in sequencing depth
    object <- logNormCounts(
        object = object, block = block, threads = threads,
        size_factors = size_factors,
        name = "multiBatchNorm"
    )

    # dimensionality reduction
    object <- runPCA(
        object = object,
        assay = "multiBatchNorm",
        threads = threads, name = "PCA",
        d = d, scale = scale, subset_row = subset_row,
        block = block, block_weight_policy = block_weight_policy,
        variable_block_weight = variable_block_weight,
        from_residuals = from_residuals, extra_work = extra_work,
        iterations = iterations, seed = seed,
        realized = realized
    )

    # run MNN for batch correction
    runMNN(
        object = object, block = block, dimred = "PCA",
        ..., threads = threads
    )
}

#' @export
#' @rdname quickMNN
quickMNN.Seurat <- function(
    object, block, size_factors = NULL,
    # `runPCA` arguments
    d = 50L, scale = FALSE, subset_row = NULL,
    block_weight_policy = NULL,
    variable_block_weight = c(0, 1000),
    from_residuals = FALSE, extra_work = 7,
    iterations = 1000, seed = NULL,
    realized = TRUE,
    # `runMNN` arguments
    ..., threads = NULL) {
    # nromalization, adjust for differences in sequencing depth
    object <- logNormCounts(
        object = object, block = block, threads = threads,
        size_factors = size_factors,
        name = "multiBatchNorm"
    )

    # dimensionality reduction
    object <- runPCA(
        object = object, block = block,
        layer = "multiBatchNorm",
        threads = threads, name = "PCA",
        d = d, scale = scale, subset_row = subset_row,
        block = block, block_weight_policy = block_weight_policy,
        variable_block_weight = variable_block_weight,
        from_residuals = from_residuals, extra_work = extra_work,
        iterations = iterations, seed = seed,
        realized = realized
    )

    # run MNN for batch correction
    runMNN(
        object = object, block = block, dimred = "PCA",
        ..., threads = threads
    )
}
