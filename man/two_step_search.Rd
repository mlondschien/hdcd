% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimizer.R
\name{two_step_search}
\alias{two_step_search}
\title{Two step search optimisation algorithm}
\usage{
two_step_search(gain_function, split_candidates, n, control)
}
\arguments{
\item{gain_function}{closure}

\item{split_candidates}{An array of indices where gain function is (possibly) evaluated.}

\item{n}{total number of observations}

\item{control}{Control parameter as returned by \link{hdcd_control}}
}
\value{
A list with entries gain, max_gain, best_split, permutation_test and pval
}
\description{
Find an (approximate) maximum of \code{gain_function} by repeatedly fitting and evaluating.
}
