`%||%` <- function(x, y) if (is.null(x)) y else x

set_threads <- function(threads,
                        arg = rlang::caller_arg(threads),
                        call = rlang::caller_call()) {
    assert_number_whole(threads,
        min = 1L, allow_null = TRUE, arg = arg, call = call
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
check_seed <- function(seed, arg = caller_arg(seed), call = caller_call()) {
    if (is.null(seed)) {
        if (is.null(old_seed())) on.exit(restore_rng(NULL))
        random_seed(1L)
    } else {
        assert_number_whole(seed, arg = arg, call = call)
        seed
    }
}

#' @importFrom rlang caller_env
set_seed <- function(seed, envir = caller_env()) {
    run <- substitute(on.exit(restore_rng(oseed)), list(oseed = old_seed()))
    eval(run, envir = envir)
    set.seed(seed)
}

old_seed <- function() {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
        NULL
    }
}

restore_rng <- function(oseed) {
    if (is.null(oseed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        }
    } else {
        assign(".Random.seed", oseed, envir = .GlobalEnv, inherits = FALSE)
    }
}

random_seed <- function(n) sample.int(1e6L, n)
