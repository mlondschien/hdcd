% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{RegrowNetwork}
\alias{RegrowNetwork}
\title{RegrowNetwork}
\usage{
RegrowNetwork(omega, n_nodes = ncol(omega) * 0.1,
  preferential_power = 1, v = 0.3)
}
\arguments{
\item{omega}{A precision matrix as created by ScaleNetwork}

\item{n_nodes}{Number of nodes to prune and regrow. Default is 0.1 of all nodes.}

\item{preferential_power}{Power coefficient alpha for weighting of degree number as alpha in prefential attachment mechanism.}

\item{v}{Constant added to the off diagonal of the precision matrix for controlling the magnitude of partial correlations.}
}
\description{
Prune graph and regrow edges according of scale-free network
}
\details{
Note that v and preferential_power agruments need to be equal to the ones which
initially created omega.
}
\examples{
omega <- ScaleNetwork(20, v = 1)
omega_hat <- RegrowNetwork(omega, v = 1)
omega
omega_hat
}
