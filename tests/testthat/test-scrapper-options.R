test_that("variance modelling supports block quantiles and optional trends", {
    d <- compat_data()
    out <- modelGeneVar(d$logcounts, block = d$block,
        block_average_policy = "quantile", block_quantile = 0.25)
    raw <- scrapper::modelGeneVariances(d$logcounts, block = d$block,
        block.average.policy = "quantile", block.quantile = 0.25,
        use.min.width = FALSE)
    expect_equal(as.matrix(out), as.matrix(raw$statistics))

    if (utils::packageVersion("scrapper") >= "1.7.1") {
        out <- modelGeneVar(d$logcounts, block = d$block, fit_trend = FALSE)
        expect_named(out, c("means", "variances"))
        expect_named(attr(out, "per_block")[[1]], c("means", "variances"))
        out <- modelGeneVar(d$logcounts, block = d$block,
            fit_trend = FALSE, block_average_policy = "none")
        expect_identical(dim(out), c(400L, 0L))
        expect_length(attr(out, "per_block"), 2L)
    }
})

test_that("PCA subsets retain the full rotation and pass convergence controls", {
    d <- compat_data()
    out <- runPCA(d$logcounts, n_dim = 5, subset = 1:40, seed = 42, warn = FALSE)
    raw <- scrapper::runPca(d$logcounts, number = 5, subset = 1:40,
        seed = 42, warn = FALSE, extra.work = 7)
    expect_equal(as.vector(out), as.vector(t(raw$components)))
    expect_equal(attr(out, "rotation"), raw$rotation)
    expect_identical(nrow(attr(out, "rotation")), 400L)

    if (utils::packageVersion("scrapper") >= "1.6.3") {
        out <- runPCA(d$logcounts, n_dim = 5, seed = 42, tolerance = 1e-7)
        raw <- scrapper::runPca(d$logcounts, number = 5, seed = 42,
            tolerance = 1e-7, extra.work = 7)
        expect_equal(as.vector(out), as.vector(t(raw$components)))
        scores <- scoreFeatureSet(d$logcounts, 1:40, seed = 42, tolerance = 1e-7)
        raw <- scrapper::scoreGeneSet(d$logcounts, set = 1:40, seed = 42,
            tolerance = 1e-7, extra.work = 7)
        expect_equal(unname(scores[1, ]), raw$scores)
    }
})

test_that("development scrapper options reach UMAP and feature scoring", {
    skip_if(utils::packageVersion("scrapper") < "1.7.1")
    d <- compat_data()
    out <- runUMAP(d$pcs, seed = c(42, 43), n_epochs = 10,
        optimization = TRUE, threads = 1)
    raw <- scrapper::runUmap(d$pcs, initialize.seed = 42, optimize.seed = 43,
        num.epochs = 10, min.dist = 0.01, num.threads.optimize = 1)
    expect_identical(out, raw)
    out <- scoreFeatureSet(d$logcounts, 1:40, seed = 42, warn = FALSE)
    raw <- scrapper::scoreGeneSet(d$logcounts, set = 1:40, seed = 42,
        warn = FALSE, extra.work = 7)
    expect_equal(unname(out[1, ]), raw$scores)
    expect_warning(invisible(scoreFeatureSet(d$logcounts, 1:40, seed = 42,
        iterations = 1, warn = TRUE)), "converge")
    expect_warning(invisible(scoreFeatureSet(d$logcounts, 1:40, seed = 42,
        iterations = 1, warn = FALSE)), NA)
})
