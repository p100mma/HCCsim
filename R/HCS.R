

#' Extract matrix of principal components of all clusters of HCD
#' 
#' @param \code{hcluster_blockPCA} A list of blockwise PCA reconstruction of clusters and subclusters as obtained and \code{subdivideNreconstruct} function.
#' @return A matrix with all prinicipal components inside each block of \code{hcluster_blockPCA} structure. Column names give tuples of clusters to which each prinicipal component belongs.
#' @examples
#' data(brca)
#' data(brca_clusters)
#' lvl1<- initial_clusterNreconstruct(X= brca, X_variances=matrixStats::colVars(brca),
#'				        clustering_vector=brca_clusters)
#' lvl2<- subclusterNreconstruct(X=brca,
#' 				 X_variances= matrixStats::colVars(brca),
#' 				 hclustering_list= lvl1$clustering_list,
#'                               hcluster_PCA_blocks=lvl1$cluster_blockPCA,
#'                               clfun2=similarity_based_hclust,
#'                               clfun2OtherArgs_constant=list(method="complete" ),
#'				 clfun2OtherArgs_ranges= list(n_group=2:7)) 
#' PC_generator_matrix( lvl2$hcluster_blockPCA) -> PCmat
#'  #note: no correlations between base clusters and subclusters.
#' heatmap( sqrt(abs(cor(PCmat))), Colv=NA, Rowv=NA)
#' @export

PC_generator_matrix<- function(hcluster_blockPCA){

all_PC<- do.call( cbind, lapply(hcluster_blockPCA, function(block) block$PC) )
colnames(all_PC)<-do.call(c, lapply( seq_along(hcluster_blockPCA), function(i) 
					rep( names(hcluster_blockPCA)[[i]],
					     hcluster_blockPCA[[i]]$k_used)
				    )
			)
all_PC
}

#' Generate samples from multivariate normal distribution with a given covariance
#' 
#' @param sigma A covariance matrix to use. Must be positive definite.
#' @param n_samples Number of samples to generate
#' @return A matrix of generated samples from a mutlivariate normal distribution with covariance \code{sigma}, each column corresponds to each variable and row to each sample.
#' @examples
#' data(brca)
#' multivariate_normal_cholesky( cov(brca[,1:300] ), 600) -> gen_data
#' dim(gen_data)
#' MSE(cov(brca[,1:300]), cov(gen_data) )
#' @export

multivariate_normal_cholesky<- function( sigma, n_samples){
U = chol(sigma)
t(U) -> L
Z<- matrix(rnorm( ncol(sigma)*n_samples), nrow= n_samples )
t( L %*% t(Z) )
}



