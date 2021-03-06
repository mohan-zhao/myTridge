% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{Gradient}
\alias{Gradient}
\title{Gradient of objective function}
\usage{
Gradient(theta, family, y, X, trex_c)
}
\arguments{
\item{theta}{unknown regression vector}

\item{family}{three family cases: gaussian, poisson and binomial}

\item{y}{where y = Xβ + ε}

\item{X}{high dimensional matrix}

\item{trex_c}{TREX parameters, 2 for gaussian case, 1 for poisson and binomial case}
}
\description{
Gradient of objective function
}
