# Standalone file: do not edit by hand
# Source: <https://github.com/Yunuuuu/standalone/blob/main/R/standalone-pkg.R>
# ----------------------------------------------------------------------
#
# ---
# repo: Yunuuuu/standalone
# file: standalone-pkg.R
# last-updated: 2025-03-04
# license: https://unlicense.org
# imports: [utils]
# ---

# This file contains various helper utilities, including common functions
# used across multiple packages I have developed. Some functions depend on
# other packages that are not listed in Imports, so use them with caution.

# ## Changelog
# 2025-03-04:
# - Add `%||%`
#
# 2025-03-03:
# - Add `rd_collect_family`
# - Add `oxford_and`
# - Add `oxford_or`
# - Add `code_quote`
# - Add `oxford_comma`
#
# 2025-02-26:
# - Add `is_installed`
# - Add `install_pkgs`
# - Add `pkg_nm`
# - Add `pkg_namespace`
#
# nocov start

`%||%` <- function(x, y) if (is.null(x)) y else x

is_installed <- local({
    cache <- new.env(parent = emptyenv())
    function(pkg, version = NULL) {
        id <- if (is.null(version)) pkg else paste(pkg, version, sep = ":")
        out <- cache[[id]]
        if (is.null(out)) {
            if (is.null(version)) {
                out <- requireNamespace(pkg, quietly = TRUE)
            } else {
                out <- requireNamespace(pkg, quietly = TRUE) &&
                    utils::packageVersion(pkg) >= version
            }
            assign(id, out, envir = cache, inherits = FALSE)
        }
        out
    }
})

install_pkgs <- function(pkgs) {
    if (is_installed("pak")) {
        getExportedValue("pak", "pkg_install")(pkgs, ask = FALSE)
    } else {
        utils::install.packages(pkgs)
    }
}

pkg_nm <- function() utils::packageName(environment())

pkg_namespace <- function() topenv(environment())

# Need `rlang` package
set_exit <- function(expr, envir = parent.frame(), after = TRUE, add = TRUE) {
    expr <- getExportedValue("rlang", "enquo")(expr)
    thunk <- as.call(list(
        getExportedValue("rlang", "new_function")(list(), expr)
    ))
    do.call(base::on.exit, list(thunk, add = add, after = after), envir = envir)
}

# utils function to collapse characters ---------------------------
oxford_and <- function(x, code = TRUE, quote = TRUE, sep = ", ") {
    oxford_comma(code_quote(x, code, quote), sep = sep, final = "and")
}

oxford_or <- function(x, code = TRUE, quote = TRUE, sep = ", ") {
    oxford_comma(code_quote(x, code, quote), sep = sep, final = "or")
}

code_quote <- function(x, code = TRUE, quote = TRUE) {
    if (quote) x <- paste0("\"", x, "\"")
    if (code) x <- paste0("`", x, "`")
    x
}

oxford_comma <- function(x, sep = ", ", final = "and") {
    n <- length(x)

    if (n < 2L) return(x) # styler: off

    head <- x[seq_len(n - 1L)]
    last <- x[n]

    head <- paste(head, collapse = sep)

    # Write a or b. But a, b, or c.
    if (n > 2L) {
        paste0(head, sep, final, " ", last)
    } else {
        paste0(head, " ", final, " ", last)
    }
}

# Need `roxygen2` package
#' @description add `@eval rd_collect_family("myfamily")` to the functions in
#' your package. This will automatically generate a section listing all
#' functions tagged with `@family myfamily`.
#' @param family A string specifying the family name.
#' @param section_title A string specifying the section title.
#' @param code_style A boolean indicating whether to apply code formatting
#' to function names.
#' @noRd
rd_collect_family <- function(family,
                              section_title = paste(family, "family"),
                              code_style = TRUE) {
    # get blocks objects from the roxygenize function
    blocks <- NULL
    pos <- sys.nframe()
    while (pos > 0L) {
        if (!is.null(call <- sys.call(-pos))) {
            fn <- eval(.subset2(call, 1L), sys.frame(-(pos + 1L)))
            env <- sys.frame(-pos)
            if (identical(fn, getExportedValue("roxygen2", "roxygenize")) &&
                exists("blocks", envir = env, inherits = FALSE)) {
                blocks <- get("blocks", envir = env, inherits = FALSE)
                break
            }
        }
        pos <- pos - 1L
    }

    # identify the blocks with family of the same tag specified in `family`
    blocks <- blocks[
        vapply(blocks, function(block) {
            getExportedValue("roxygen2", "block_has_tags")(block, "family") &&
                identical(
                    getExportedValue("roxygen2", "block_get_tag_value")(
                        block, "family"
                    ),
                    family
                )
        }, logical(1L), USE.NAMES = FALSE)
    ]
    if (length(blocks) == 0L) return(character()) # styler: off

    # extracted the function name
    funs <- vapply(blocks, function(block) {
        as.character(.subset2(block$call, 2L))
    }, character(1L), USE.NAMES = FALSE)
    if (code_style) {
        items <- sprintf("\\code{\\link[=%s]{%s()}}", funs, funs)
    } else {
        items <- sprintf("\\link[=%s]{%s()}", funs, funs)
    }
    c(
        sprintf("@section %s:", section_title),
        "\\itemize{",
        sprintf("  \\item %s", items),
        "}"
    )
}

# nocov end
