% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BRCA_GE-data.R
\docType{data}
\name{brca}
\alias{brca}
\title{BRCA gene expression dataset}
\format{
A large numeric matrix where columns are genes and rows are samples.
}
\source{
\url{https://ega-archive.org/studies/EGAS00000000083}
}
\usage{
data(brca)
}
\description{
A subset of 500 random samples and 2500 random genes from BRCA dataset retrieved from European Genome-Phenome Archive, accesion ID: EGAS00000000083.
}
\examples{
data(brca)
dim(brca)
}
\references{
Pereira B, Chin   SF, Rueda   OM  et al.   The somatic mutation profiles of 2,433 breast cancers refine their genomic and transcriptomic landscapes. Nat Commun  2016;7:11479
}
\seealso{
\link{brca_clusters} for precomputed clustering of genes in this data
}
\keyword{datasets}
