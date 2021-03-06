% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun.r
\name{dcmp}
\alias{dcmp}
\alias{pcmp}
\alias{rcmp}
\title{Conway-Maxwell-Poisson Distribution}
\usage{
dcmp(x, lognu, logmu, maxiter = 50)

pcmp(q, lognu, logmu, maxiter = 50)

rcmp(logmu, lognu, maxiter = 50)
}
\arguments{
\item{x}{vector of (non-negative integer) quantiles.}

\item{lognu}{the shape parameter using Guikema and Goffelt's (2008) parametrization.}

\item{logmu}{the mean using Guikema and Goffelt's (2008) parametrization.}

\item{maxiter}{integer specifying number of term to keep in calculating the normalizing constant. Defaults to 50.}

\item{q}{vector of quantiles.}
}
\description{
Density, distribution function and random generaton for the Conway-Maxwell-Poisson distribution with parameters \code{logmu} and \code{lognu}.
}
\section{Functions}{
\itemize{
\item \code{dcmp}: Density function for the CMP distribution.

\item \code{pcmp}: Distribution function for the CMP distribution.

\item \code{rcmp}: Generates random numbers from the CMP distribution. This is a wrapper for the k-aggregated CMP random number generator.
}}

\references{
[1] Guikema, S. D., & Goffelt, J. P. (2008). A flexible count data regression model for risk analysis. Risk Analysis: An International Journal, 28(1), 213-223.
}
