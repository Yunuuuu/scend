`%||%` <- function(x, y) if (is.null(x)) y else x

set_threads <- function(threads, arg = rlang::caller_arg(threads),
                        call = rlang::caller_call()) {
    assert_number_whole(threads,
        min = 1, allow_null = TRUE, arg = arg, call = call
    )
    if (is.null(threads)) {
        parallel::detectCores()
    } else {
        as.integer(threads)
    }
}

is_scalar_numeric <- function(x) length(x) == 1L && is.numeric(x)

quickdf <- function(x) {
    class(x) <- "data.frame"
    attr(x, "row.names") <- .set_row_names(length(.subset2(x, 1L)))
    x
}

#' @importFrom rlang caller_arg caller_call
arg_match <- function(x, values, default = .subset(values, 1L),
                      arg = caller_arg(x), call = caller_call()) {
    if (is.null(x)) {
        default
    } else {
        rlang::arg_match0(
            arg = x, values = values,
            arg_nm = arg, error_call = call
        )
    }
}
