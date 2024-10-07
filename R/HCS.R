

#' Extract matrix of principal components of all clusters of HCD
#' 
#' @param \code{hcluster_blockPCA} A list of blockwise PCA reconstruction of clusters and subclusters as obtained and \code{subdivideNreconstruct} function.
#' @return A matrix with all prinicipal components inside each block of \code{hcluster_blockPCA} structure. Column names give tuples of clusters to which each prinicipal component belongs.
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
#' @export

multivariate_normal_cholesky<- function( sigma, n_samples){
U = chol(sigma)
t(U) -> L
Z<- matrix(rnorm( ncol(sigma)*n_samples), nrow= n_samples )
t( L %*% t(Z) )
}



