#' Clusters of genes of BRCA dataset
#'
#' Clusters from a run of Markov Clustering algorithm [2] on BRCA dataset, as was done in publication [1]. Clustering was performed on the whole dataset whereas cluster labels of this vector correspond to random subset of 2500 genes selected for \code{brca} data. For details go to [1].
#'
#' @docType data
#'
#' @usage data(brca_clusters)
#'
#' @format A vector of integers encoding cluster labels of each gene. Zeroes denote clusterless genes.
#'
#' @keywords cluster_labels
#' 
#' @references [1]  Piotr Stomma, Witold R Rudnicki, HCS—hierarchical algorithm for simulation of omics datasets,  Bioinformatics, Volume 40, Issue Supplement_2, September 2024, Pages ii98–ii104, \url{https://doi.org/10.1093/bioinformatics/btae392}
#' @references [2] van Dongen, Stijn, Graph clustering via a discrete uncoupling process, Siam Journal on Matrix Analysis and Applications 30-1, p121-141, 2008
#'
#' @seealso [brca] for the data corresponding to the clusters
#' @examples 
#' data(brca_clusters)
#' length(brca_clusters)
#' table(brca_clusters)
"brca_clusters"
