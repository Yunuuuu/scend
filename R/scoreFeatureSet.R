#' Score feature set activity for each cell
#'
#' @inheritParams runPCA
#' @export
scoreFeatureSet <- function(object, ...) {
    UseMethod("scoreFeatureSet")
}

#' @param feature_sets A list of integer, logical or character vector specifying
#' the features that belong to the set. You can also directly input an atomic
#' vector if there is only one feature set.
#' @inheritParams runPCA
#' @inheritParams scrapper::scoreGeneSet
#' @inherit scrapper::scoreGeneSet description details
#' @param warn Whether to warn if IRLBA does not converge. Requires scrapper
#' 1.7.1 or later when set explicitly.
#' @return A matrix of feature scores. The `weights` attribute is a list with
#' one data frame per feature set, containing `row` (the feature's row index in
#' `object`) and `weight` columns.
#' @seealso [`scoreGeneSet`][scrapper::scoreGeneSet]
#' @rdname scoreFeatureSet
#' @export
scoreFeatureSet.default <- function(object, feature_sets, ...,
                                    rank = 1, scale = FALSE,
                                    block = NULL, block_weight_policy = NULL,
                                    variable_block_weight = c(0, 1000),
                                    extra_work = 7, iterations = 1000,
                                    seed = NULL, realized = TRUE,
                                    threads = NULL, tolerance = NULL,
                                    warn = TRUE) {
    rlang::check_dots_empty()
    threads <- set_threads(threads)
    if (!is.list(feature_sets)) feature_sets <- list(feature_sets)
    if (!length(feature_sets)) {
        return(structure(
            matrix(numeric(),
                nrow = 0L, ncol = ncol(object),
                dimnames = list(NULL, colnames(object))
            ),
            weights = list()
        ))
    }
    seed <- check_seed(seed, length(feature_sets))
    irlba_args <- list()
    if (!is.null(tolerance)) irlba_args$tolerance <- tolerance
    if (utils::packageVersion("scrapper") >= "1.7.1" || !missing(warn)) {
        irlba_args$warn <- warn
    }
    out <- .mapply(function(features, seed) {
        .scoreFeatureSet(
            object = object,
            features = features,
            rank = rank, scale = scale,
            block = block, block_weight_policy = block_weight_policy,
            variable_block_weight = variable_block_weight,
            extra_work = extra_work,
            iterations = iterations,
            irlba_args = irlba_args,
            seed = seed,
            realized = realized, threads = threads
        )
    }, list(features = feature_sets, seed = seed), NULL)
    names(out) <- names(feature_sets)
    weights <- lapply(out, function(x) {
        as.data.frame(attr(x, "weights", exact = TRUE))
    })
    out <- rlang::inject(base::rbind(!!!out))
    rownames(out) <- names(feature_sets) %||% seq_along(feature_sets)
    colnames(out) <- colnames(object)
    structure(out, weights = weights)
}

#' @export
scoreFeatureSet.HDF5Matrix <- function(object, ..., threads = NULL) {
    threads <- set_threads(threads)
    if (threads > 1L) {
        cli::cli_warn(c(
            "Cannot use multiple threads for {.cls HDF5Matrix}",
            i = "Will use {.code threads = 1} instead"
        ))
        threads <- 1L
    }
    NextMethod()
}

#' @param name A string of the assay name for the scores.
#' @export
#' @rdname scoreFeatureSet
scoreFeatureSet.SummarizedExperiment <- function(object, ...,
                                                 name = "scores",
                                                 assay = "logcounts") {
    mat <- .get_mat_from_sce(object, assay)
    scores <- list(scoreFeatureSet(object = mat, ...))
    names(scores) <- name
    SummarizedExperiment::SummarizedExperiment(
        assays = scores,
        colData = SummarizedExperiment::colData(object)
    )
}

#' @export
#' @rdname scoreFeatureSet
scoreFeatureSet.Seurat <- function(object, ..., assay = NULL, layer = "data") {
    mat <- .get_mat_from_seurat(object, assay, layer)
    scoreFeatureSet(object = mat, ...)
}

.scoreFeatureSet <- function(object, features, rank = 1, scale = FALSE,
                             block = NULL, block_weight_policy = NULL,
                             variable_block_weight = c(0, 1000),
                             extra_work = 7,
                             iterations = 1000, seed = 1234L,
                             realized = TRUE, threads = NULL,
                             irlba_args = list()) {
    out <- rlang::inject(scrapper::scoreGeneSet(
        object,
        set = features,
        rank = rank,
        scale = scale,
        block = block,
        block.weight.policy = block_weight_policy,
        variable.block.weight = variable_block_weight,
        extra.work = extra_work,
        iterations = iterations,
        seed = seed,
        realized = realized,
        num.threads = threads,
        !!!irlba_args
    ))
    structure(.subset2(out, "scores"), weights = .subset2(out, "weights"))
}
