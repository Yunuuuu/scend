`%||%` <- function(x, y) if (is.null(x)) y else x

set_threads <- function(threads,
                        arg = rlang::caller_arg(threads),
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

#' @importFrom rlang caller_arg caller_call
check_seed <- function(seed, len = 1L,
                       arg = caller_arg(seed), call = caller_call()) {
    if (is.null(seed)) {
        if (is.null(old_seed())) on.exit(restore_rng(NULL))
        seed <- random_seed(1L)
    } else if (len == 1L) {
        assert_number_whole(seed, arg = arg, call = call)
    }
    if (len > 1L) seed <- check_seeds(seed, len, arg = arg, call = call)
    seed
}

check_seeds <- function(seed, len, arg = caller_arg(seed),
                        call = caller_call()) {
    if (length(seed) == 1L) {
        set_seed(seed)
        random_seed(len)
    } else if (length(seed) != len) {
        cli::cli_abort(
            "{.arg {arg}} must be a scalar or of length {len}",
            call = call
        )
    } else {
        seed[seq_len(len)]
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
