#' Find subclusters
#'
#' @inheritParams runPCA
#' @export
subCluster <- function(object, ...) UseMethod("subCluster")

#' @inheritParams scoreMarkers
#' @param restricted Character vector containing the subset of groups in
#' `groups` to be subclustered. By default, all unique groups in `groups` are
#' used for subclustering, but this can be restricted to specific groups of
#' interest to save compute time.
#' @param ... Additional arguments passed on to `clusterFn`.
#' @param clusterFn A function returns a vector of cluster assignments for each
#' cell in that object.
#' @param format A string to be passed to [`sprintf`], specifying how
#' the subclusters should be named with respect to the parent level in
#' `groups` and the level returned by `clusterFn`.
#' @param new_level A string of `"insert"`, `"fill-end"`, `"start"`, and `"end"`
#' indicates how to deal with the new levels. `"insert"` means the new levels
#' will be added into the old group level, `"fill-start"` and `"fill-end"` means
#' the first of the new levels will be insert into the old group level position
#' and the redundant new levels will be added in the start or end respectively.
#' `"start"` and `"end"` means the new levels will be added in the start or end,
#' and the old group level will be removed directly. Default: `"fill-end"`.
#' @return An integer vector with cluster assignments for each cell. Additional
#' attributes return by `clusterFn` will be added for each subcluster run.
#' @export
#' @rdname subCluster
subCluster.default <- function(object, groups, restricted = NULL, ...,
                               clusterFn = clusterSNNGraph, format = "%s.%s",
                               new_level = "fill-end") {
    # Seurat::FindSubCluster()
    policy <- rlang::arg_match0(
        new_level,
        c("insert", "fill-end", "fill-start", "start", "end")
    )
    by_groups <- split(seq_along(groups), groups)
    if (!is.null(restricted)) {
        restricted <- as.character(restricted)
        missing <- setdiff(restricted, names(by_groups))
        if (length(missing)) {
            cli::cli_abort(c(
                "Invalid value found in {.arg restricted}",
                i = "Cannot find {.val {missing}} group{?s}"
            ))
        }
        by_groups <- .subset(by_groups, restricted)
    }
    # Do subclustering for each group ------------------------
    all_groups <- as.character(groups)
    all_levels <- levels(groups) %||% sort(unique(all_groups))
    all_attributes <- NULL
    for (group in names(by_groups)) {
        index <- .subset2(by_groups, group)
        sub_membership <- clusterFn(object = object[, index], ...)
        all_groups[index] <- sprintf(format, group, sub_membership)
        # deal with the new-generated levels ---------------------
        new_levels <- sprintf(format, group, levels(sub_membership))
        all_levels <- order_levels(policy, group, all_levels, new_levels)
        # save the attributes returned by `clusterFn` ------------
        if (!is.null(sub_attributes <- attributes(sub_membership))) {
            names(sub_attributes) <- sprintf(
                "subCluster%s_%s", group, names(sub_attributes)
            )
            all_attributes <- c(all_attributes, sub_attributes)
        }
    }
    out <- factor(all_groups, levels = all_levels)
    for (attr_nm in names(all_attributes)) {
        attr(out, attr_nm) <- all_attributes[[attr_nm]]
    }
    out
}

order_levels <- function(policy, current, old, new) {
    i <- which(old == current)
    if (policy == "insert") {
        append(old[-i], new, i - 1L)
    } else if (policy == "fill-end") {
        old[i] <- new[1L]
        c(old, new[-1L])
    } else if (policy == "fill-start") {
        old[i] <- new[1L]
        c(new[-1L], old)
    } else if (policy == "end") {
        c(old[-i], new)
    } else if (policy == "start") {
        c(new, old[-i])
    }
}
