% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics.R
\name{SE}
\alias{SE}
\title{Compute squared error between two vectors}
\usage{
SE(x, y)
}
\arguments{
\item{x}{One of vectors between which squared error is to be computed.}

\item{y}{One of vectors between which squared error is to be computed.}
}
\value{
A vector of squared errors between \code{x,y}.
}
\description{
Compute squared error between two vectors
}
\examples{
SE( c(1,2,3), c(3,2,1) )
}
