% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph_construction.R
\name{similarity_matrix}
\alias{similarity_matrix}
\title{Compute similarity matrix from input correlation coefficients}
\usage{
similarity_matrix(
  corData,
  addLoops = TRUE,
  power = 2,
  level = 0.05,
  return.lower.tri = FALSE
)
}
\arguments{
\item{corData}{A list storing lower triangular halves of correlation and p-value matrices, see output of \code{fastPearsonData} or \code{rcorrData} for how it should be structured.}

\item{addLoops}{Logical, if \code{TRUE}, diagonal of the similarity matrix will assume values of \code{1} nad \code{0} otherwise.}

\item{power}{A power to which raise the absolute value of the correlation coefficients.}

\item{level}{A significance level of adjusted p-values over which similarities are set to zero.}

\item{return.lower.tri}{A logical, if \code{TRUE} then function returns lower triangular half of the similarity matrix. Otherwise, full (symmetric) matrix of similarities is returned.}
}
\value{
Depending on the value of \code{return.lower.tri}, either a square, symmetric matrix of nonnegative similarities between objects ( \code{S[i,j]= abs(cor(X_i,X_j))^power} if \code{cor(X_i,X_j)} is significant at \code{level} after p-value adjustment), or a lower triangular half of that matrix (\code{S[lower.tri(S)]})
}
\description{
Function takes compressed data about symmetric correlation matrix and builds a \code{WGCNA} style similarity matrix out of it, raising elements to \code{power} and zeroing out unsignificant correlations.
}
\examples{
data(brca)
brca_corData= fastPearsonData(brca)
S= similarity_matrix( brca_corData)
ncol(S)==ncol(brca)
}
