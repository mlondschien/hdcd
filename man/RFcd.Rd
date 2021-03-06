% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hdcd.R
\name{RFcd}
\alias{RFcd}
\title{Random Forest change point detection}
\usage{
RFcd(x, delta = 0.1, control = hdcd_control())
}
\arguments{
\item{x}{A matrix with observations in rows}

\item{delta}{Minimal relative segment length, defaults to 0.1}

\item{control}{Control parameter as returned by \link{hdcd_control}}
}
\description{
Find (non-parametric) breaks in the distribution of a time series using Random Forest Classifiers
}
