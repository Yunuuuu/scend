% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runUMAP.R
\name{runUMAP}
\alias{runUMAP}
\alias{runUMAP.default}
\alias{runUMAP.SingleCellExperiment}
\alias{runUMAP.Seurat}
\title{Compute the uniform manifold approximation and projection with \code{scrapper}}
\usage{
runUMAP(object, ...)

\method{runUMAP}{default}(
  object,
  n_dim = 2L,
  n_neighbors = 15L,
  n_epochs = -1L,
  min_dist = 0.01,
  ...,
  optimization = FALSE,
  seed = NULL,
  bnparam = AnnoyParam(),
  threads = NULL
)

\method{runUMAP}{SingleCellExperiment}(
  object,
  ...,
  dimred = "PCA",
  n_dimred = NULL,
  assay = NULL,
  name = "UMAP"
)

\method{runUMAP}{Seurat}(
  object,
  ...,
  dimred = "PCA",
  n_dimred = NULL,
  assay = NULL,
  layer = NULL,
  name = "UMAP"
)
}
\arguments{
\item{object}{A matrix-like object containing the counts (non-negative
integers). Rows are features and columns are cells.}

\item{...}{\itemize{
\item \code{default} method: not used.
\item \code{SingleCellExperiment} and \code{Seurat} method: additional arguments passed to
\code{default} method.
}}

\item{n_dim}{Integer scalar specifying the number of dimensions of the output
embedding.}

\item{n_neighbors}{Integer scalar specifying the number of neighbors to use
in the UMAP algorithm.}

\item{n_epochs}{Integer scalar specifying the number of epochs to perform.
If set to \code{-1}, an appropriate number of epochs is chosen based on
\code{ncol(object)}.}

\item{min_dist}{Numeric scalar specifying the minimum distance between
points.}

\item{optimization}{Logical scalar specifying whether to parallelize the
optimization step.}

\item{seed}{Integer scalar specifying the seed to use for generating the initial coordinates.}

\item{bnparam}{A \linkS4class{BiocNeighborParam} object specifying the
algorithm to use.}

\item{threads}{Integer scalar specifying the number of threads to use. If
\code{NULL}, a single thread will be used.}

\item{dimred}{String or integer scalar specifying the existing dimensionality
reduction results to use.}

\item{n_dimred}{Integer scalar or vector specifying the dimensions to use if
dimred is specified.}

\item{assay}{Integer scalar or string indicating which assay of x
contains the expression values.}

\item{name}{String specifying the name to be used to store the result in the
\code{\link[SingleCellExperiment:reducedDims]{reducedDims}} or
\code{\link[SeuratObject:Seurat-class]{reductions}} of the output.}

\item{layer}{Name of the layer to get from the assay data.}
}
\value{
\itemize{
\item \code{default} method: A numeric matrix where rows are cells and columns are
the two dimensions of the embedding.
\item \code{SingleCellExperiment} method: embedding was added into
\code{\link[SingleCellExperiment:reducedDims]{reducedDims}} named as \code{name}.
\item \code{Seurat} method: embedding was added into
\code{\link[SeuratObject:Seurat-class]{reductions}} named as \code{name}.
}
}
\description{
Compute the uniform manifold approximation and projection with \code{scrapper}
}
\seealso{
\code{\link[scrapper:runUmap]{runUmap}}
}
