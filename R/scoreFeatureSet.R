#' Score feature set activity for each cell
#'
#' @inheritParams runPCA
#' @export
scoreFeatureSet <- function(object, ...) {
    rlang::check_dots_used()
    UseMethod("scoreFeatureSet")
}

#' @param feature_sets A list of integer, logical or character vector specifying
#' the features that belong to the set. You can also directly input an atomic
#' vector if there is only one feature set.
#' @inheritParams runPCA
#' @inheritParams scrapper::scoreGeneSet
#' @inherit scrapper::scoreGeneSet description details
#' @return A matrix of feature scores.
#' @seealso [scoreGeneSet][scrapper::scoreGeneSet]
#' @rdname scoreFeatureSet
#' @export
scoreFeatureSet.default <- function(object, feature_sets, ...,
                                    rank = 1, scale = FALSE,
                                    block = NULL, block_weight_policy = NULL,
                                    variable_block_weight = c(0, 1000),
                                    extra_work = 7,
                                    iterations = 1000, seed = NULL,
                                    realized = TRUE, threads = NULL) {
    threads <- set_threads(threads)
    seed <- check_seed(seed)
    set_seed(seed)
    if (is.list(feature_sets)) {
        out <- lapply(feature_sets, function(features) {
            .scoreFeatureSet(
                object = object,
                features = features,
                rank = rank, scale = scale,
                block = block, block_weight_policy = block_weight_policy,
                variable_block_weight = variable_block_weight,
                extra_work = extra_work,
                iterations = iterations,
                seed = random_seed(1L),
                realized = realized, threads = threads
            )
        })
        weights <- lapply(out, attr, which = "weights", exact = TRUE)
        weights <- do.call(base::cbind, weights)
        out <- do.call(base::rbind, out)
        colnames(out) <- colnames(object)
        structure(out, weights = weights)
    } else {
        out <- .scoreFeatureSet(
            object = object,
            features = feature_sets,
            rank = rank, scale = scale,
            block = block, block_weight_policy = block_weight_policy,
            variable_block_weight = variable_block_weight,
            extra_work = extra_work,
            iterations = iterations,
            seed = random_seed(1L),
            realized = realized, threads = threads
        )
        dim(out) <- c(1L, length(out))
        colnames(out) <- colnames(object)
        out
    }
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
