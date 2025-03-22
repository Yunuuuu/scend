
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Efficient End-to-End Analysis for Single-Cell RNA-seq data (scend)

<!-- badges: start -->
<!-- badges: end -->

## Installation

You can install the development version of scend from
[GitHub](https://github.com/) with:

``` r
if (!requireNamespace("pak")) {
    install.packages("pak",
        repos = sprintf(
            "https://r-lib.github.io/p/pak/devel/%s/%s/%s",
            .Platform$pkgType, R.Version()$os, R.Version()$arch
        )
    )
}
pak::pkg_install("Yunuuuu/scend")
```

## Introduction

[scrapper](https://github.com/libscran/scrapper) is a package provides
methods for an end-to-end analysis of a single-cell RNA-sequencing
(scRNA-seq) analysis, starting from the count matrix and finishes with
clusters, markers, and various embeddings (i.e., t-SNE and UMAP). It’s
pretty fast and memory-efficient. The package serves as a bridging tool
connecting the `scrapper` package with either `Seurat` or
`SingleCellExperiment` objects.

## Available functions

| functions            | Description                                             |
|----------------------|---------------------------------------------------------|
| `librarySizeFactors` | Compute library size factors                            |
| `logNormCounts`      | Log-transformed normalized expression                   |
| `modelGeneVar`       | Model per-gene variances in expression                  |
| `runPCA`             | Principal component analysis                            |
| `runMNN`             | Fast mutual nearest neighbors correction                |
| `subSample`          | Subsample cells based on their neighbors                |
| `runTSNE`            | t-stochastic neighbor embedding                         |
| `runUMAP`            | uniform manifold approximation and projection           |
| `runDiffusionMap`    | diffusion map                                           |
| `clusterSNNGraph`    | Graph-based clustering                                  |
| `subCluster`         | Find subclusters                                        |
| `scoreMarkers`       | Score marker genes                                      |
| `scoreFeatureSet`    | Score feature set activity for each cell                |
| `infercnv`           | Run infercnv workflow                                   |
| `scenic`             | Single-Cell Regulatory Network Inference and Clustering |
|                      |                                                         |
