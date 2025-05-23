% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scoreMarkers.R
\name{scoreMarkers}
\alias{scoreMarkers}
\alias{scoreMarkers.SingleCellExperiment}
\alias{scoreMarkers.Seurat}
\alias{scoreMarkers.default}
\title{Score marker genes}
\usage{
scoreMarkers(object, ...)

\method{scoreMarkers}{SingleCellExperiment}(object, ..., assay = "logcounts")

\method{scoreMarkers}{Seurat}(object, ..., assay = NULL, layer = "data")

\method{scoreMarkers}{default}(
  object,
  groups,
  ...,
  block = NULL,
  block_weight_policy = NULL,
  variable_block_weight = c(0, 1000),
  compute_auc = TRUE,
  threshold = 0L,
  all_pairwise = FALSE,
  threads = NULL
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

\item{assay}{Integer scalar or string indicating which assay of x
contains the expression values.}

\item{layer}{Name of the layer to get from the assay data.}

\item{groups}{A vector specifying the group assignment for each cell in
\code{object}.}

\item{block}{Factor specifying the block of origin (e.g., batch, sample) for
each cell in \code{x}. Alternatively \code{NULL} if all cells are from the same block.}

\item{block_weight_policy}{A string of \code{"variable"}, \code{"equal"}, or \code{"none"}
specifying the policy to use for weighting different blocks when computing
the average for each statistic Only used if \code{block} is not \code{NULL}.}

\item{variable_block_weight}{Numeric vector of length 2, specifying the
parameters for variable block weighting. The first and second values are used
as the lower and upper bounds, respectively, for the variable weight
calculation. Only used if \code{block} is not \code{NULL} and \code{block_weight_policy = "variable"}.}

\item{compute_auc}{Logical scalar indicating whether to compute the \code{AUC}.
Setting this to \code{FALSE} can improve speed and memory efficiency.}

\item{threshold}{Non-negative numeric scalar specifying the minimum threshold
on the differences in means (i.e., the log-fold change, if \code{object} contains
log-expression values). This is incorporated into the effect sizes for
\verb{Cohen's d} and the \code{AUC}.}

\item{all_pairwise}{Logical scalar indicating whether to report the full
effects for every pairwise comparison.}

\item{threads}{Integer scalar specifying the number of threads to use. If
\code{NULL}, a single thread will be used.}
}
\value{
A list of data frame of marker statistics.
Each data frame corresponds to a group in \code{groups} and contains:
\itemize{
\item \code{mean}: the mean expression across all cells in the current group.
\item \code{detected}: proportion of cells with detectable expression in the current
group.
\item \code{cohens.d}: the Cohen's d statistics across all pairwise comparisons
involving the current group. This includes the \code{min}, \code{mean},
\code{median}, \code{max} and \code{rank} (\code{min.rank}).
\item \code{auc}: the AUC statistics across all pairwise comparisons
involving the current group. This includes the \code{min}, \code{mean},
\code{median}, \code{max} and \code{rank} (\code{min.rank}).
\item \code{delta.mean}: the difference in the mean expression compared to other
groups. This includes the \code{min}, \code{mean}, \code{median}, \code{max} and
\code{rank} (\code{min.rank}).
\item \code{delta.detected}: the difference in the detected proportions compared to
other groups. This includes the \code{min}, \code{mean}, \code{median},
\code{max} and \code{rank} (\code{min.rank}).
}

If \code{all_pairwise=TRUE}, this list will also contain \code{pairwise}, a list of
lists of data frames. Each data frame contains the statistics for the
pairwise comparison between groups, e.g., \code{pairwise$A$B} contains the
statistics for \verb{A versus B} where large effects correspond to upregulation in
\code{A}.
}
\description{
Score marker genes
}
\seealso{
\itemize{
\item \link[scrapper:scoreMarkers]{scrapper::scoreMarkers}
\item \link[scrapper:summarizeEffects]{summarizeEffects}
}
}
