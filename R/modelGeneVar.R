#' @inherit scrapper::modelGeneVariances title description
#' @param object A matrix-like object containing the counts (non-negative
#' integers). Rows are features and columns are cells.
#' @param ...
#'  - `default` method: not used.
#'  - `SingleCellExperiment` and `Seurat` method: additional arguments passed to
#'    `default` methods.
#' @return A data frame of class `scend_modelGeneVar` with the columns `means`,
#' `variances`, `fitted` `residuals`, each of which is a numeric vector
#' containing the statistic of the same name across all genes.
#'
#' If `block` is supplied, each of the column vectors described above contains
#' the average across all blocks. The data frame will also contain `per_block`
#' attribute, a list of data frames containing the equivalent statistics for
#' each block.
#' @export
modelGeneVar <- function(object, ...) {
    UseMethod("modelGeneVar")
}

#' @param assay Integer scalar or string indicating which assay of x
#' contains the expression values.
#' @export
#' @rdname modelGeneVar
modelGeneVar.SingleCellExperiment <- function(object, ...,
                                              assay = "logcounts") {
    mat <- .get_mat_from_sce(object, assay)
    modelGeneVar(object = mat, ...)
}

#' @param layer Name of the layer to get from the assay data.
#' @export
#' @rdname modelGeneVar
modelGeneVar.Seurat <- function(object, ..., assay = NULL, layer = "data") {
    mat <- .get_mat_from_seurat(object, assay, layer)
    modelGeneVar(object = mat, ...)
}

#' @param block Factor specifying the block of origin (e.g., batch, sample) for
#' each cell in `x`. Alternatively `NULL` if all cells are from the same block.
#' @param block_weight_policy A string of `"variable"`, `"equal"`, or `"none"`
#' specifying the policy to use for weighting different blocks when computing
#' the average for each statistic Only used if `block` is not `NULL`.
#' @param variable_block_weight Numeric vector of length 2, specifying the
#' parameters for variable block weighting. The first and second values are used
#' as the lower and upper bounds, respectively, for the variable weight
#' calculation. Only used if `block` is not `NULL` and `block_weight_policy =
#' "variable"`.
#' @param mean_filter Logical scalar indicating whether to filter on the means
#' before trend fitting.
#' @param min_mean Numeric scalar specifying the minimum mean of genes to use in
#' trend fitting.  Only used if `mean_filter=TRUE`.
#' @param transform Logical scalar indicating whether a quarter-root
#' transformation should be applied before trend fitting.
#' @param span Numeric scalar specifying the span of the LOWESS smoother.
#' Ignored if `use_min_width=TRUE`.
#' @param use_min_width Logical scalar indicating whether a minimum width
#' constraint should be applied to the LOWESS smoother.  Useful to avoid
#' overfitting in high-density intervals.
#' @param min_width Minimum width of the window to use when
#' `use_min_width=TRUE`.
#' @param min_window_count Minimum number of observations in each window. Only
#' used if `use_min_width=TRUE`.
#' @param threads Integer scalar specifying the number of threads to use. If
#' `NULL`, all detected threads will be used. See
#' [`detectCores`][parallel::detectCores].
#' @export
#' @rdname modelGeneVar
modelGeneVar.default <- function(object,
                                 block = NULL, block_weight_policy = NULL,
                                 ...,
                                 variable_block_weight = c(0, 1000),
                                 mean_filter = TRUE, min_mean = 0.1,
                                 transform = TRUE, span = 0.3,
                                 use_min_width = FALSE, min_width = 1,
                                 min_window_count = 200, threads = NULL) {
    rlang::check_dots_empty()
    threads <- set_threads(threads)
    if (!is.null(block)) {
        block_weight_policy <- rlang::arg_match0(
            block_weight_policy, c("variable", "equal", "none")
        )
    }
    stats <- scrapper::modelGeneVariances(
        x = object,
        block = block, block.weight.policy = block_weight_policy,
        variable.block.weight = variable_block_weight,
        mean.filter = mean_filter,
        min.mean = min_mean, transform = transform,
        span = span, use.min.width = use_min_width,
        min.width = min_width, min.window.count = min_window_count,
        num.threads = threads
    )
    ans <- .subset2(stats, "statistics")
    attr(ans, "per_block") <- .subset2(stats, "per.block")
    class(ans) <- c("scend_modelGeneVar", class(ans))
    ans
}
