% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim.R
\name{sim}
\alias{sim}
\title{Monte Carlo simulation for the filter estimation}
\usage{
sim(N, lambda, mu, delta, alpha, kappa, filt = c("cutoff", "landw"),
  m0 = 0, alphaLoss = -1)
}
\arguments{
\item{N}{Integer number of simulation runs.}

\item{lambda}{Numeric vector of decreasing, strictly positive entries of
the diagonal design matrix.}

\item{mu}{Numeric vector valued input signal.}

\item{delta}{Numeric noise level.}

\item{alpha}{Numeric smoothing index for the residuals.}

\item{kappa}{Numeric critical stopping value.}

\item{filt}{Character string designating the filter to be used. filt should
be one of "cutoff" or "landw".}

\item{m0}{Integer cut-off index for the two-step procudure.}

\item{alphaLoss}{Numeric smoothing parameter for the loss. \code{alpha} =
-1 gives the strong loss. \code{alpha} = 0 gives the weak loss.}
}
\value{
Returns a data frame containing \code{N} observations of the
 residual stopping rule and loss.
}
\description{
Performs \code{N} simulation runs of the filter estimation procedure for a
given design matrix \code{lambda}, a given signal \code{mu} and a given
noise level \code{delta} and evaluates the loss and oracle quantities.
}
