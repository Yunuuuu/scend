#' @importFrom rlang caller_arg caller_call
check_seed <- function(seed, len = 1L,
                       arg = caller_arg(seed), call = caller_call()) {
    if (is.null(seed)) {
        restore_rng_hook()
        seed <- random_seed(len)
    } else if (len == 1L) {
        assert_number_whole(seed, arg = arg, call = call)
    } else if (length(seed) == 1L) { # we need multiple seeds
        restore_rng_hook()
        set.seed(seed)
        seed <- random_seed(len)
    } else if (length(seed) != len) {
        cli::cli_abort(
            "{.arg {arg}} must be a single number or of length {len}",
            call = call
        )
    } else {
        seed <- seed[seq_len(len)]
    }
    seed
}

#' @importFrom rlang caller_env
set_seed <- function(seed, envir = caller_env()) {
    restore_rng_hook(envir = envir)
    set.seed(seed)
}

old_rng <- function() {
    # styler: off
    if (exists(".Random.seed", envir = .GlobalEnv,
               mode = "integer",  inherits = FALSE)) {
        # styler: on
        get(".Random.seed",
            envir = .GlobalEnv,
            mode = "integer", inherits = FALSE
        )
    } else {
        NULL
    }
}

restore_rng <- function(orng) {
    if (is.null(orng)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        }
    } else {
        assign(".Random.seed", orng, envir = .GlobalEnv, inherits = FALSE)
    }
}

restore_rng_hook <- function(orng = old_rng(), envir = parent.frame()) {
    code <- substitute(on.exit(restore_rng(orng)), list(orng = orng))
    eval(code, envir = envir)
}

random_seed <- function(n) sample.int(1e6L, n)

seed <- function(x = NULL, rng_kind = NULL, rng_normal_kind = NULL) {
    structure(
        list(seed = x, rng_kind = rng_kind, rng_normal_kind = rng_normal_kind),
        class = "scend_seed"
    )
}
