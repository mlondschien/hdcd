% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gain_functions.R
\name{glasso_gain_function}
\alias{glasso_gain_function}
\title{Closure generating function to calculate gains when splitting and learning a glasso model}
\usage{
glasso_gain_function(x, start, end, lambda, control = hdcd_control())
}
\arguments{
\item{x}{A matrix with observations in rows}

\item{start}{The start of the investigated interval}

\item{end}{The end of the investigated interval}

\item{lambda}{Tuning parameter passed to the method.}

\item{control}{an object of type \code{hdcd_control} returned by \link{hdcd_control}}
}
\value{
A closure with parameters \code{x}, \code{start} and \code{end}, that when evaluated
will itself return a closure with parameter \code{split_point}. This calculates the gain when
splitting the segment \code{(start, end]} of \code{x} at \code{split_point}. If the closure is
additionally supplied with \code{evaluate_all = TRUE}, an array of length \code{nrow(x)} is
returned with differences of loglikelihoods for each observation in \code{(start, end]} when split at
\code{split_point}.
Glasso Gain Function

A closure with parameters \code{x}, \code{start} and \code{end}, that when evaluated
will itself return a closure with arguments \code{split_point} and \code{split_candidates}.
Evaluated this function returns the gain when splitting the segment \code{(start, end]} of \code{x}
 at \code{split_point}.
}
\description{
Closure generating function to calculate gains when splitting and learning a glasso model
}
