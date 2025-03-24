#' @importFrom rlang caller_arg caller_call
check_seed <- function(seed, len = 1L,
                       arg = caller_arg(seed), call = caller_call()) {
    if (is.null(seed)) {
    } else if (len == 1L) {
        assert_number_whole(seed, arg = arg, call = call)
    } else if (length(seed) == 1L) { # we need multiple seeds
        restore_seed_hook()
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

#' @param seed An integer nubmer.
#' @importFrom rlang caller_env
#' @noRd
set_seed <- function(seed, envir = caller_env()) {
    if (!is.null(seed)) {
        restore_seed_hook(envir = envir)
        set.seed(seed)
    }
}

has_seed <- function() {
    exists(".Random.seed",
        envir = globalenv(),
        mode = "integer", inherits = FALSE
    )
}

get_seed <- function() {
    if (!has_seed()) {
        return(NULL)
    }
    list(
        rng = get(".Random.seed",
            envir = globalenv(),
            mode = "integer", inherits = FALSE
        ),
        rng_kind = RNGkind()
    )
}

rm_seed <- function() {
    if (!has_seed()) {
        return(NULL)
    }
    set.seed(seed = NULL)
    rm(".Random.seed", envir = globalenv())
}

restore_seed <- function(seed) {
    if (is.null(seed)) {
        rm_seed()
    } else {
        restore_rng_kind(seed$rng_kind)
        if (is.null(seed$rng)) {
            if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
                rm(".Random.seed", envir = globalenv(), inherits = FALSE)
            }
        } else {
            set.seed(seed$rng)
        }
    }
}

#' @importFrom rlang caller_env
restore_seed_hook <- function(after = TRUE, seed = get_seed(),
                              envir = caller_env()) {
    set_exit(restore_seed(!!seed), envir = envir, after = after)
}

restore_rng_kind <- function(rng_kind) {
    RNGkind(.subset2(rng_kind, 1L), normal.kind = .subset2(rng_kind, 2L))
    if (identical(sample_kind <- .subset2(rng_kind, 3L), "Rounding")) {
        suppressWarnings(RNGkind(sample.kind = sample_kind))
    } else {
        RNGkind(sample.kind = sample_kind)
    }
}

random_seed <- function(n) sample.int(1e6L, n)

seed <- function(x = NULL, rng_kind = NULL, rng_normal_kind = NULL) {
    structure(
        list(rng = x, rng_kind = rng_kind, rng_normal_kind = rng_normal_kind),
        class = "scend_seed"
    )
}
