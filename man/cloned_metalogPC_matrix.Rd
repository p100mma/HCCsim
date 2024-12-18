% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HCS.R
\name{cloned_metalogPC_matrix}
\alias{cloned_metalogPC_matrix}
\title{Get a set of cloned PCs transformed to target metalog distributions based on input normal PCs}
\usage{
cloned_metalogPC_matrix(
  PChidden_normal,
  target_metalogs,
  target_terms,
  n_samples = nrow(PChidden_normal)
)
}
\arguments{
\item{target_metalogs}{A list of \code{metalog} objects specifying target marginal distributions for generated variables. Must be of length equal to number of columns of \code{PC_hidden_normal}.}

\item{target_terms}{A vector of numbers of terms to use in applying each metalog transformation of the same length as \code{target_metalogs}.}

\item{n_samples}{Number of samples to generate, defaults to the number of rows of \code{reference_PC_gen_matrix}.}

\item{PC_hidden_normal}{A matrix of which covariance is computed and replicated in output data. Column names of that matrix get copied to the output one.}
}
\value{
A matrix of same dimension as \code{PC_hidden_normal} containing cloned version of input matrix which columns are uncorrelated with the input and are transformed to target metalog distributions.
}
\description{
Get a set of cloned PCs transformed to target metalog distributions based on input normal PCs
}
\details{
Function will generate set of variables having same covariance as columns of \code{PChidden_normal}. Then,
each such variable gets transformed to metalog distribution \link{1} given by \code{target_metalogs[[i]]}
(a \code{metalog} object).
}
\examples{
data(brca)
data(brca_clusters)
lvl1<- initial_clusterNreconstruct(X= brca, X_variances=matrixStats::colVars(brca),
			        clustering_vector=brca_clusters)
PC_generator_matrix( lvl1$cluster_blockPCA) -> PCmat
meta_list<- list()
for (j in 1:ncol(PCmat))  rmetalog::metalog(PCmat[,j], term_limit=5, step_len=.01) -> meta_list[[j]]
HN_PCmat<- hidden_normalPC_matrix(PCmat, meta_list, rep(5, length(meta_list) ) )
cloned_metaPC<- cloned_metalogPC_matrix(HN_PCmat, meta_list, rep(5, length(meta_list)) )
 #note: no correlations between cloned and reference PCs
heatmap( abs(cor(cbind(PCmat, cloned_metaPC))), Colv=NA, Rowv=NA)
}
\references{
\link{1} Keelin T.  The metalog distributions. Dec Anal  2016;13:243–77.
}
\seealso{
\code{\link[=hidden_normalPC_matrix]{hidden_normalPC_matrix()}}, \code{\link[=PC_generator_matrix]{PC_generator_matrix()}} ,\code{\link[=subclusterNreconstruct]{subclusterNreconstruct()}} ,\code{\link[=initial_clusterNreconstruct]{initial_clusterNreconstruct()}} for how to generate input data for this function, \code{\link[rmetalog:metalog]{rmetalog::metalog()}} for metalog distribution
}
