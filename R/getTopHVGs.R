#' @inherit scrapper::chooseHighlyVariableGenes title description
#' @inheritParams rlang::args_dots_used
#' @param stats Numeric vector of variances (or a related statistic) across all
#' genes. Typically the `residuals` from [`modelGeneVar`] are used here.
#' @inheritParams scrapper::chooseHighlyVariableGenes
#' @param keep_ties Logical scalar indicating whether to keep tied values of
#' `stats`, even if `top` may be exceeded.
#' @return A character of highly variable genes.
#' @export
getTopHVGs <- function(stats, ...) {
    rlang::check_dots_used()
    UseMethod("getTopHVGs")
}

#' @export
#' @rdname getTopHVGs
getTopHVGs.scend_modelGeneVar <- function(stats, top = 2000,
                                          larger = TRUE,
                                          keep_ties = TRUE, ...) {
    assert_number_whole(top)
    assert_bool(larger)
    assert_bool(keep_ties)
    scrapper::chooseHighlyVariableGenes(
        .subset2(stats, "residuals"),
        top = top, larger = larger, keep.ties = keep_ties
    )
}
