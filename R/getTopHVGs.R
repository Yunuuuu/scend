#' @inherit scrapper::chooseHighlyVariableGenes title description
#' @inheritParams rlang::args_dots_used
#' @param n Integer specifying the number of top genes to retain.
#' @param prop Numeric scalar specifying the proportion of genes to report as
#' HVGs.
#' @param stats A data frame returned by [`modelGeneVar`]. Its `residuals`
#' column is used to select highly variable genes.
#' @inheritParams scrapper::chooseHighlyVariableGenes
#' @param keep_ties Logical scalar indicating whether to keep tied values of
#' `stats`, even if `n` may be exceeded.
#' @return An integer vector of row indices for the highly variable genes.
#' @export
getTopHVGs <- function(stats, ...) {
    UseMethod("getTopHVGs")
}

#' @export
#' @rdname getTopHVGs
getTopHVGs.scend_modelGeneVar <- function(stats, n = NULL, prop = NULL,
                                          larger = TRUE, keep_ties = TRUE,
                                          ...) {
    rlang::check_dots_empty()
    assert_number_whole(n, min = 0, allow_null = TRUE)
    assert_number_decimal(prop, min = 0, max = 1, allow_null = TRUE)
    assert_bool(larger)
    assert_bool(keep_ties)
    if (!is.null(n) && !is.null(prop)) {
        n <- max(n, round(prop * nrow(stats)))
    } else if (!is.null(n)) {} else if (!is.null(prop)) {
        n <- round(prop * nrow(stats))
    } else {
        n <- nrow(stats)
    }
    scrapper::chooseHighlyVariableGenes(
        .subset2(stats, "residuals"),
        top = n,
        larger = larger,
        keep.ties = keep_ties
    )
}
