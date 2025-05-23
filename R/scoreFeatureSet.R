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
#' @return A matrix of feature scores.
#' @seealso [`scoreGeneSet`][scrapper::scoreGeneSet]
#' @rdname scoreFeatureSet
#' @export
scoreFeatureSet.default <- function(object, feature_sets, ...,
                                    rank = 1, scale = FALSE,
                                    block = NULL, block_weight_policy = NULL,
                                    variable_block_weight = c(0, 1000),
                                    extra_work = 7, iterations = 1000,
                                    seed = NULL, realized = TRUE,
                                    threads = NULL) {
    rlang::check_dots_empty()
    threads <- set_threads(threads)
    if (!is.list(feature_sets)) feature_sets <- list(feature_sets)
    seed <- check_seed(seed, length(feature_sets))
    out <- .mapply(function(features, seed) {
        .scoreFeatureSet(
            object = object,
            features = features,
            rank = rank, scale = scale,
            block = block, block_weight_policy = block_weight_policy,
            variable_block_weight = variable_block_weight,
            extra_work = extra_work,
            iterations = iterations,
            seed = seed,
            realized = realized, threads = threads
        )
    }, list(features = feature_sets, seed = seed), NULL)
    weights <- lapply(out, attr, which = "weights", exact = TRUE)
    weights <- do.call(base::cbind, weights)
    out <- do.call(base::rbind, out)
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
                             realized = TRUE, threads = NULL) {
    out <- scrapper::scoreGeneSet(
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
        num.threads = threads
    )
    structure(.subset2(out, "scores"), weights = .subset2(out, "weights"))
}
