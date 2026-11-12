test_that("default seeds produce usable and reproducible embeddings", {
    d <- compat_data()
    set.seed(123)
    pca <- runPCA(d$logcounts, n_dim = 5)
    set.seed(123)
    expect_identical(runPCA(d$logcounts, n_dim = 5), pca)
    expect_identical(dim(pca), c(60L, 5L))
    expect_warning(runPCA(d$logcounts, n_dim = 5, block = d$block,
        from_residuals = TRUE, seed = 123), NA)

    set.seed(123)
    tsne <- runTSNE(d$pcs, perplexity = 5, max_iter = 10)
    set.seed(123)
    expect_identical(runTSNE(d$pcs, perplexity = 5, max_iter = 10), tsne)
    expect_identical(dim(tsne), c(60L, 2L))
    expected <- suppressWarnings(scrapper::runTsne(d$pcs, perplexity = 5,
        num.neighbors = 12, max.iterations = 10, max.depth = 20, seed = 42))
    actual <- suppressWarnings(runTSNE(d$pcs, perplexity = 5,
        n_neighbors = 12, max_iter = 10, seed = 42))
    expect_identical(actual, expected)
})

test_that("UMAP maps seeds and automatic epochs to supported upstream options", {
    x <- compat_data()$pcs
    expect_warning(automatic <- runUMAP(x, seed = 42), NA)
    expect_identical(dim(automatic), c(60L, 2L))
    expect_identical(runUMAP(x, n_epochs = -1, seed = 42), automatic)
    expect_identical(runUMAP(x, seed = 42), automatic)
    expect_warning(short <- runUMAP(x, n_epochs = 10,
        seed = c(42, 43), optimization = TRUE, threads = 1), NA)
    expect_identical(dim(short), c(60L, 2L))
    expect_identical(runUMAP(x, n_epochs = 10,
        seed = c(42, 43), optimization = FALSE, threads = 1), short)
    expect_identical(dim(runUMAP(x, n_epochs = 10)), c(60L, 2L))
})

test_that("MNN exposes effective parameters", {
    d <- compat_data()
    expect_warning(out <- runMNN(d$pcs, d$block,
        n_steps = 2, merge_policy = "input"), NA)
    raw <- scrapper::correctMnn(d$pcs, d$block, num.steps = 2, merge.policy = "input")
    expect_equal(as.vector(out), as.vector(t(raw$corrected)))
    expect_identical(rownames(out), colnames(d$pcs))
})

test_that("quick workflows work with default seeds and explicit threads", {
    d <- compat_data()
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = d$counts))
    expect_warning(corrected <- quickMNN(sce, d$block, n_dim = 5, threads = 1), NA)
    expect_identical(dim(SingleCellExperiment::reducedDim(corrected, "corrected")), c(60L, 5L))
    clusters <- quickCluster(sce, n_dim = 5, n_hvgs = 40,
        n_neighbors = 5, method = "louvain", threads = 1)
    expect_length(clusters, ncol(sce))
    expect_false(anyNA(clusters))
})
