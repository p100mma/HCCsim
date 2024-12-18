% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_utils.R
\name{clustering_vector2list}
\alias{clustering_vector2list}
\title{Transform clustering_vector to clustering_list}
\usage{
clustering_vector2list(clustering_vector)
}
\arguments{
\item{clustering_vector}{An integer vector for which \code{clustering_vector[i]} gives a label of the group to which node \code{i} belongs. Its attribute \code{'params'} is copied onto the resulting \code{clustering_list}. Does not have to be an object of class \code{HCCSim_clustering_vector} but function assumes it sticks to the labeling convention outlined for that kind of representation.}
}
\value{
A list with named elements with one element per proper cluster, where each element is integer vector of indexes of nodes. An object of class \code{HCCSim_clustering_list}.
}
\description{
Function converts a vector of numeric labels to the \code{clustering_list} representation of the clusters. See fifunctions \code{HCCSim_clustering_vector} and \code{HCCSim_clustering_list}.
}
\examples{
data(brca_clusters)
clustering_vector2list( HCCSim_clustering_vector(brca_clusters, params =list(algorithm='mcl', inflation=2)) )
}
\seealso{
\code{\link[=HCCSim_clustering_vector]{HCCSim_clustering_vector()}}, \code{\link[=HCCSim_clustering_list]{HCCSim_clustering_list()}} and \code{\link[=clustering_list2vector]{clustering_list2vector()}}
}
