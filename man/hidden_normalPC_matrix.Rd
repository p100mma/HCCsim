% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HCS.R
\name{hidden_normalPC_matrix}
\alias{hidden_normalPC_matrix}
\title{Generate matrix of PCs transformed to normal distribution}
\usage{
hidden_normalPC_matrix(reference_PC_gen_matrix, target_metalogs, target_terms)
}
\arguments{
\item{reference_PC_gen_matrix}{Matrix where columns are PCs of clusters which will be transformed to normal distributions.}

\item{target_metalogs}{A list of objects of class \code{metalog} of length equal to number of columns of \code{reference_PC_gen_matrix}.}

\item{target_terms}{A vector of numbers of terms to use in applying each metalog transformation of the same length as \code{target_metalogs}.}
}
\value{
A matrix of the same dimensions as \code{reference_PC_gen_matrix}, containing transformed variables from input matrix, with same column names.
}
\description{
Generate matrix of PCs transformed to normal distribution
}
\details{
First, each \code{reference_PC_gen_matrix[,i]} gets transformed by applying CDF of metalog \link{1} distribution
stored in \code{target_metalogs[i]} ( \code{metalog} object) using \code{target_terms[i]} terms.
Then, variable number \code{i} is pulled under standard normal distribution by applying inverse of its
cumulative distribution function.
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
#display the effect of the transform
par(mfrow=c(1,2))
hist(PCmat[,1])
hist(HN_PCmat[,1])
}
\references{
\link{1} Keelin T.  The metalog distributions. Dec Anal  2016;13:243–77.
}
\seealso{
\code{\link[=PC_generator_matrix]{PC_generator_matrix()}} ,\code{\link[=subclusterNreconstruct]{subclusterNreconstruct()}} ,\code{\link[=initial_clusterNreconstruct]{initial_clusterNreconstruct()}} for how to generate input data for this function, \code{\link[rmetalog:metalog]{rmetalog::metalog()}} for metalog distribution
}
