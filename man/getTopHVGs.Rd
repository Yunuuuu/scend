% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getTopHVGs.R
\name{getTopHVGs}
\alias{getTopHVGs}
\alias{getTopHVGs.scend_modelGeneVar}
\title{Choose highly variable genes}
\usage{
getTopHVGs(stats, ...)

\method{getTopHVGs}{scend_modelGeneVar}(stats, n = NULL, prop = NULL, larger = TRUE, keep_ties = TRUE, ...)
}
\arguments{
\item{stats}{Numeric vector of variances (or a related statistic) across all
genes. Typically the \code{residuals} from \code{\link{modelGeneVar}} are used here.}

\item{...}{Arguments passed to methods.}

\item{n}{Integer specifying the number of top genes to retain.}

\item{prop}{Numeric scalar specifying the proportion of genes to report as
HVGs.}

\item{larger}{Logical scalar indicating whether larger values of \code{stats} correspond to more variable genes.
If \code{TRUE}, HVGs are defined as those with the largest values of \code{stats}.}

\item{keep_ties}{Logical scalar indicating whether to keep tied values of
\code{stats}, even if \code{n} may be exceeded.}
}
\value{
A character of highly variable genes.
}
\description{
Choose highly variable genes (HVGs) based on a variance-related statistic.
}
