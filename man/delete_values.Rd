% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{delete_values}
\alias{delete_values}
\title{Delete values from a design matrix}
\usage{
delete_values(x, m, missingness = "mcar")
}
\arguments{
\item{x}{Design matrix}

\item{m}{proportion of values to be deleted}

\item{missingness}{One of \code{'mcar'}, \code{'nmar'} or \code{'blockwise'}. For \code{mcar} the values are
deleted uniformly at random. If \code{nmar} is supplied, the probability of deletion of a given value is set 
to be proportional to its absolute value. If \code{blockwise} is supplied, repeatedly values \code{l ~ Pois(p / 20)}
and \code{k ~ Exp(n/8)} are drawn and a block of size \code{l x k} is deleted.}
}
\description{
\code{delete_values} returns the matrix \code{x} with a proportion of \code{m} values deleted.
}
