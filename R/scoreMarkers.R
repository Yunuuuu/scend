#' Score marker genes
#' @export
scoreMarkers <- function(object, ...) {
    UseMethod("scoreMarkers")
}

#' @export
#' @rdname scoreMarkers
scoreMarkers.SingleCellExperiment <- function(object, ..., assay = "logcounts") {
    mat <- .get_mat_from_sce(object, assay)
    scoreMarkers(object = mat, ...)
}

#' @export
#' @rdname scoreMarkers
scoreMarkers.Seurat <- function(object, ..., assay = NULL, layer = "data") {
    mat <- .get_mat_from_seurat(object, assay, layer)
    scoreMarkers(object = mat, ...)
}

#' @inheritParams modelGeneVar
#' @param groups A vector specifying the group assignment for each cell in
#' `object`.
#' @param compute_auc Logical scalar indicating whether to compute the `AUC`.
#' Setting this to `FALSE` can improve speed and memory efficiency.
#' @param threshold Non-negative numeric scalar specifying the minimum threshold
#' on the differences in means (i.e., the log-fold change, if `object` contains
#' log-expression values). This is incorporated into the effect sizes for
#' `Cohen's d` and the `AUC`.
#' @param all_pairwise Logical scalar indicating whether to report the full
#' effects for every pairwise comparison.
#' @inheritParams logNormCounts
#' @inherit scrapper::scoreMarkers details return
#' @return A list of data frame of marker statistics.
#' Each data frame corresponds to a group in `groups` and contains:
#'
#' - `mean`: the mean expression across all cells in the current group.
#' - `detected`: proportion of cells with detectable expression in the current
#'             group.
#' - `cohens.d`: the Cohen's d statistics across all pairwise comparisons
#'              involving the current group. This includes the `min`, `mean`,
#'              `median`, `max` and `rank` (`min.rank`).
#' - `auc`: the AUC statistics across all pairwise comparisons
#'              involving the current group. This includes the `min`, `mean`,
#'              `median`, `max` and `rank` (`min.rank`).
#' - `delta.mean`: the difference in the mean expression compared to other
#'                 groups. This includes the `min`, `mean`, `median`, `max` and
#'                `rank` (`min.rank`).
#' - `delta.detected`: the difference in the detected proportions compared to
#'                     other groups. This includes the `min`, `mean`, `median`,
#'                    `max` and `rank` (`min.rank`).
#'
#' If `all_pairwise=TRUE`, this list will also contain `pairwise`, a list of
#' lists of data frames. Each data frame contains the statistics for the
#' pairwise comparison between groups, e.g., `pairwise$A$B` contains the
#' statistics for `A versus B` where large effects correspond to upregulation in
#' `A`.
#'
#' @seealso
#' - [scrapper::scoreMarkers]
#' - [summarizeEffects][scrapper::summarizeEffects]
#' @rdname scoreMarkers
#' @export
scoreMarkers.default <- function(object, groups,
                                 ...,
                                 block = NULL, block_weight_policy = NULL,
                                 variable_block_weight = c(0, 1000),
                                 compute_auc = TRUE,
                                 threshold = 0L, all_pairwise = FALSE,
                                 threads = NULL) {
    rlang::check_dots_empty()
    threads <- set_threads(threads)
    scores <- scrapper::scoreMarkers(
        x = object,
        groups = groups,
        block = block,
        block.weight.policy = block_weight_policy,
        variable.block.weight = variable_block_weight,
        compute.auc = compute_auc,
        threshold = threshold,
        all.pairwise = all_pairwise,
        num.threads = threads
    )
    # re-style the scores to the same with `scran.chan::scoreMarkers.chan()`
    stats <- c("cohens.d", "auc", "delta.mean", "delta.detected")
    lvls <- levels(groups) %||% sort(unique(groups))
    names(lvls) <- lvls
    if (all_pairwise) {
        pairwise <- lapply(lvls, function(e) { # experimental group
            lapply(lvls, function(c) { # control group
                ans <- vector("list", length(stats))
                names(ans) <- stats
                ans
            })
        })
        for (i in stats) {
            array <- .subset2(scores, i)
            for (e in lvls) {
                for (c in lvls) {
                    pairwise[[e]][[c]][[i]] <- array[e, c, , drop = TRUE]
                }
            }
            scores[[i]] <- scrapper::summarizeEffects(array, threads)
        }
    }
    out <- lapply(lvls, function(group) {
        cur_stats <- vector("list", 2L + length(stats))
        names(cur_stats) <- c("mean", "detected", stats)
        cur_stats$mean <- .subset2(scores, "mean")[
            , group,
            drop = TRUE
        ]
        cur_stats$detected <- .subset2(scores, "detected")[
            , group,
            drop = TRUE
        ]
        for (i in stats) {
            stat <- .subset2(.subset2(scores, i), group)
            cur_stats[[i]] <- rlang::set_names(stat, function(nm) {
                sub("min\\.rank$", "rank", nm)
            })
        }
        do.call(cbind, cur_stats)
    })
    if (all_pairwise) {
        out$pairwise <- lapply(pairwise, function(data) {
            lapply(data, function(l) {
                class(l) <- "data.frame"
                attr(l, "row.names") <- names(.subset2(l, 1L))
                l
            })
        })
    }
    out
}

#' @export
scoreMarkers.HDF5Matrix <- function(object, ..., threads = NULL) {
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
