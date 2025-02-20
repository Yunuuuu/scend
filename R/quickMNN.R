#' Quick fastMNN with `scrapper`
#'
#' @inheritDotParams runMNN
#' @inheritParams runPCA
#' @param lognorm_args A list of additional arguments passed on to
#' [`logNormCounts()`].
#' @inheritParams logNormCounts
#' @seealso
#' - [`logNormCounts`]
#' - [`runPCA`]
#' - [`runMNN`]
#' @export
quickMNN <- function(object, ...) UseMethod("quickMNN")

#' @export
#' @rdname quickMNN
quickMNN.SingleCellExperiment <- function(object, block,
                                          # `logNormCounts` arguments
                                          size_factors = NULL, mode = NULL,
                                          # `runPCA` arguments
                                          n_dim = 50L, scale = FALSE,
                                          subset_row = NULL,
                                          block_weight_policy = NULL,
                                          variable_block_weight = c(0, 1000),
                                          from_residuals = FALSE,
                                          extra_work = 7,
                                          iterations = 1000, seed = NULL,
                                          realized = TRUE,
                                          # `runMNN` arguments
                                          ..., threads = NULL,
                                          # additional arguments
                                          lognorm_args = list()) {
    # nromalization, adjust for differences in sequencing depth --------
    # This can be done before variance-modelling or after.
    # - https://github.com/LTLA/batchelor/blob/master/R/quickCorrect.R; in the
    #   source code of `quickCorrect`, variance-modelling is performed after
    #   multiBatchNorm.
    # All following threads suggested run variance-modelling in the original
    # counts.
    # - https://bioconductor.riken.jp/packages/3.8/workflows/vignettes/simpleSingleCell/inst/doc/work-5-mnn.html
    # - https://support.bioconductor.org/p/9145895/#9145905
    # - http://bioconductor.org/books/3.16/OSCA.multisample/integrating-datasets.html#slower-setup
    # So I prefered to use multiBatchNorm only for the batch corrected logcounts
    # (Don't touch the original logcounts). If the computer has a large memory,
    # it's better to use `dgCMatrix` since it's much faster and accurater.
    #
    # - multiBatchNorm need size factor, if it's NULL, the internal will
    #   calculate it with `librarySizeFactors` function
    # - For normalizetion, I prefer never to use the `subset.row` since library
    #   size is the total sum of counts across all genes for each cell
    lognorm_args$block <- block
    lognorm_args$name <- lognorm_args$name %||% "multiBatchNorm"
    lognorm_args$threads <- lognorm_args$threads %||% threads
    lognorm_args$size_factors <- size_factors %||% lognorm_args$size_factors
    lognorm_args$mode <- mode %||% lognorm_args$size_factors
    lognorm_args$object <- NULL
    object <- rlang::inject(logNormCounts(object = object, !!!lognorm_args))

    # dimensionality reduction
    object <- runPCA(
        object = object,
        assay = lognorm_args$name,
        threads = threads, name = "PCA",
        n_dim = n_dim, scale = scale, subset_row = subset_row,
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
quickMNN.Seurat <- function(object, block,
                            # `logNormCounts` arguments
                            size_factors = NULL, mode = NULL,
                            # `runPCA` arguments
                            n_dim = 50L, scale = FALSE, subset_row = NULL,
                            block_weight_policy = NULL,
                            variable_block_weight = c(0, 1000),
                            from_residuals = FALSE, extra_work = 7,
                            iterations = 1000, seed = NULL,
                            realized = TRUE,
                            # `runMNN` arguments
                            ..., threads = NULL,
                            # additional arguments
                            lognorm_args = list()) {
    # nromalization, adjust for differences in sequencing depth --------
    # This can be done before variance-modelling or after.
    # - https://github.com/LTLA/batchelor/blob/master/R/quickCorrect.R; in the
    #   source code of `quickCorrect`, variance-modelling is performed after
    #   multiBatchNorm.
    # All following threads suggested run variance-modelling in the original
    # counts.
    # - https://bioconductor.riken.jp/packages/3.8/workflows/vignettes/simpleSingleCell/inst/doc/work-5-mnn.html
    # - https://support.bioconductor.org/p/9145895/#9145905
    # - http://bioconductor.org/books/3.16/OSCA.multisample/integrating-datasets.html#slower-setup
    # So I prefered to use multiBatchNorm only for the batch corrected logcounts
    # (Don't touch the original logcounts). If the computer has a large memory,
    # it's better to use `dgCMatrix` since it's much faster and accurater.
    #
    # - multiBatchNorm need size factor, if it's NULL, the internal will
    #   calculate it with `librarySizeFactors` function
    # - For normalizetion, I prefer never to use the `subset.row` since library
    #   size is the total sum of counts across all genes for each cell
    lognorm_args$block <- block
    lognorm_args$name <- lognorm_args$name %||% "multiBatchNorm"
    lognorm_args$threads <- lognorm_args$threads %||% threads
    lognorm_args$size_factors <- size_factors %||% lognorm_args$size_factors
    lognorm_args$mode <- mode %||% lognorm_args$size_factors
    lognorm_args$object <- NULL
    object <- rlang::inject(logNormCounts(object = object, !!!lognorm_args))

    # dimensionality reduction
    object <- runPCA(
        object = object,
        layer = lognorm_args$name,
        threads = threads, name = "PCA",
        n_dim = n_dim, scale = scale, subset_row = subset_row,
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
