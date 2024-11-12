`%||%` <- function(x, y) if (is.null(x)) y else x

set_threads <- function(threads,
                        arg = rlang::caller_arg(threads),
                        call = rlang::caller_call()) {
    assert_number_whole(threads, min = 1L, allow_null = TRUE)
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
