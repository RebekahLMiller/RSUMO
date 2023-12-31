% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score_sequences.R
\name{score_sequences}
\alias{score_sequences}
\title{Score sequences}
\usage{
score_sequences(
  fimo_results,
  method = c("sum_pvalue", "max_pvalue", "sum_score", "max_score")
)
}
\arguments{
\item{fimo_results}{\code{data.frame}. The table of FIMO results to use. See
\code{\link[=preprocess_fimo_results]{preprocess_fimo_results()}} for expected columns.}

\item{method}{\code{character(1)}. One of "sum_pvalue", "max_pvalue", "sum_score",
or "max_score". The scoring method to use. (Default: "sum_pvalue")}
}
\value{
A \code{data.frame} with columns \code{sequence_name}, \code{motif_id}, and \code{score}.
}
\description{
Calculate a single score for each sequence for each individual motif.
}
\examples{
print("")
}
