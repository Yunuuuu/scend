test_that("QC preserves cell names and subset proportions", {
    x <- compat_data()$counts
    qc <- perCellQCMetrics(x, subsets = list(mito = seq_len(20)))
    expect_s3_class(qc, "data.frame")
    expect_identical(rownames(qc), colnames(x))
    expect_identical(names(qc), c("sum", "detected", "subsets_mito"))
    expect_equal(qc$sum, unname(colSums(x)))
    expect_equal(qc$detected, unname(colSums(x > 0)))
    expect_equal(qc$subsets_mito, unname(colSums(x[1:20, ]) / colSums(x)))
    expect_equal(ncol(perCellQCMetrics(x)), 2)

    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = x))
    sce <- addPerCellQCMetrics(sce, subsets = list(mito = seq_len(20)))
    expect_equal(as.data.frame(SummarizedExperiment::colData(sce)), qc)
})

test_that("variance modelling preserves data frames and HVG proportions", {
    d <- compat_data()
    fit <- modelGeneVar(d$logcounts, block = d$block)
    expect_s3_class(fit, "data.frame")
    expect_identical(rownames(fit), rownames(d$logcounts))
    expect_equal(nrow(fit), 400)
    per_block <- attr(fit, "per_block")
    expect_length(per_block, 2)
    for (b in per_block) {
        expect_s3_class(b, "data.frame")
        expect_identical(rownames(b), rownames(fit))
    }
    expected <- scrapper::chooseHighlyVariableGenes(fit$residuals, top = 40)
    expect_identical(getTopHVGs(fit, prop = 0.1), expected)
    expect_identical(getTopHVGs(fit, n = 20, prop = 0.1), expected)
    expect_null(attr(modelGeneVar(d$logcounts), "per_block"))
})

test_that("marker summaries agree with upstream for numeric group labels", {
    d <- compat_data()
    for (groups in list(rep(0:2, each = 20), rep(c(2, 10, 20), each = 20),
        factor(d$groups, levels = c("C", "A", "B", "unused")))) {
        raw <- scrapper::scoreMarkers(d$logcounts, groups, compute.auc = FALSE)
        out <- scoreMarkers(d$logcounts, groups, compute_auc = FALSE)
        expect_identical(names(out), colnames(raw$mean))
        for (group in names(out)) {
            expect_s3_class(out[[group]], "data.frame")
            expect_identical(rownames(out[[group]]), rownames(d$logcounts))
            expect_equal(out[[group]]$mean, unname(raw$mean[, group]))
            expect_equal(out[[group]]$cohens.d.mean, raw$cohens.d[[group]]$mean)
            expect_false(any(startsWith(names(out[[group]]), "auc")))
        }
    }
})

test_that("pairwise marker tables retain gene identities without AUC", {
    d <- compat_data()
    for (rn in list(rownames(d$logcounts), NULL)) {
        rownames(d$logcounts) <- rn
        raw <- scrapper::scoreMarkers(d$logcounts, d$groups,
            compute.auc = FALSE, all.pairwise = TRUE)
        out <- scoreMarkers(d$logcounts, d$groups,
            compute_auc = FALSE, all_pairwise = TRUE, threads = 2)
        pair <- out$pairwise$A$B
        expect_s3_class(pair, "data.frame")
        expect_equal(nrow(pair), 400)
        expect_equal(pair$cohens.d, raw$cohens.d["A", "B", ])
        expect_identical(names(pair), c("cohens.d", "delta.mean", "delta.detected"))
        if (!is.null(rn)) expect_identical(rownames(pair), rn)
        expected <- scrapper::summarizeEffects(raw$cohens.d, num.threads = 2)
        expect_equal(out$A$cohens.d.min, expected$A$min)
    }
    with_auc <- scoreMarkers(d$logcounts, d$groups, all_pairwise = TRUE)
    expect_true("auc" %in% names(with_auc$pairwise$A$B))
})

test_that("feature scores support unequal sets and preserve weight identities", {
    d <- compat_data()
    sets <- list(first = 1:20, second = 21:50)
    scores <- scoreFeatureSet(d$logcounts, sets, seed = c(11, 12))
    expect_identical(dim(scores), c(2L, 60L))
    expect_identical(rownames(scores), names(sets))
    weights <- attr(scores, "weights")
    expect_identical(names(weights), names(sets))
    for (i in seq_along(sets)) {
        expected <- scrapper::scoreGeneSet(d$logcounts, set = sets[[i]], seed = 10 + i,
            extra.work = 7)
        expect_equal(unname(scores[i, ]), expected$scores)
        expect_equal(weights[[i]], as.data.frame(expected$weights))
        expect_equal(weights[[i]]$row, sets[[i]])
    }
    expect_identical(dim(scoreFeatureSet(d$logcounts, sets)), c(2L, 60L))
    expect_identical(dim(scoreFeatureSet(d$logcounts, 1:20)), c(1L, 60L))
    expect_identical(dim(scoreFeatureSet(d$logcounts, list())), c(0L, 60L))
})
