% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{compare_change_points}
\alias{compare_change_points}
\title{Rand type performance indices}
\usage{
compare_change_points(cpts_a, cpts_b, n)
}
\arguments{
\item{cpts_a}{A sequence of changepoints.}

\item{cpts_b}{A sequence of changepoints.}

\item{n}{Total size of dataset from which both changepoint estimates originate.}
}
\value{
Returns a vector of the index values.
}
\description{
Calculate rand type performance indices for two sets of changepoints. Typically one
of them will be the oracle estimate. See clues package for more details.
}
