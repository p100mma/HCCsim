

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

#' Get set generator PCs correlated like reference ones
#'
#' @param reference_PC_gen_matrix A matrix of PCs of clusters and subclusters which correlations are to be replicated, such as one prouced by \code{PC_generator_matrix}
#' @return A matrix of the same dimensions as input one, containing synthetic PCs from normal distribution with identical covariance as reference ones but independent of original ones.
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
#' cloned_PC<- cloned_normalPC_matrix(PCmat)
#'  #note: no correlations between cloned and reference PCs
#' heatmap( abs(cor(cbind(PCmat, cloned_PC))), Colv=NA, Rowv=NA)
#' @export


cloned_normalPC_matrix<- function(reference_PC_gen_matrix) {
	PCsig<- cov(reference_PC_gen_matrix)
	PCnorm<-multivariate_normal_cholesky(sigma=PCsig, n_samples=nrow(reference_PC_gen_matrix) ) 
	colnames(PCnorm)<- colnames(reference_PC_gen_matrix)
	PCnorm
}

#' @importFrom rmetalog pmetalog

hidden_normalPC<- function(ref_PC, target_metalog, term) {

qnorm(  pmetalog( m=target_metalog, q=ref_PC, term=term)  , 0,1)

}

#' Generate matrix of PCs transformed to normal distribution
#'
#' @details
#' First, each \code{reference_PC_gen_matrix[,i]} gets transformed by applying CDF of metalog distribution
#' stored in \code{target_metalogs[i]} ( \code{metalog} object) using \code{target_terms[i]} terms.
#' Then, variable number \code{i} is pulled under standard normal distribution by applying inverse of its 
#' cumulative distribution function.
#'
#' @param reference_PC_gen_matrix Matrix where columns are PCs of clusters which will be transformed to normal distributions. 
#' @param target_metalogs A list of objects of class \code{metalog} of length equal to number of columns of \code{reference_PC_gen_matrix}.
#' @param target_terms A vector of numbers of terms to use in applying each metalog transformation of the same length as \code{target_metalogs}.
#' @return A matrix of the same dimensions as \code{reference_PC_gen_matrix}, containing transformed variables from input matrix.
#' @examples
#' data(brca)
#' data(brca_clusters)
#' lvl1<- initial_clusterNreconstruct(X= brca, X_variances=matrixStats::colVars(brca),
#'				        clustering_vector=brca_clusters)
#' PC_generator_matrix( lvl1$cluster_blockPCA) -> PCmat
#' meta_list<- list()
#' for (j in 1:ncol(PCmat))  rmetalog::metalog(PCmat[,j], term_limit=5, step_len=.01) -> meta_list[[j]]
#' HN_PCmat<- hidden_normalPC_matrix(PCmat, meta_list, rep(5, length(meta_list) ) )
#' #display the effect of the transform
#' par(mfrow=c(1,2))
#' hist(PCmat[,1])
#' hist(HN_PCmat[,1])
#' @export

hidden_normalPC_matrix<- function(reference_PC_gen_matrix, target_metalogs, target_terms) {
stopifnot( ncol(reference_PC_gen_matrix) == length(target_metalogs) )
stopifnot( length(target_metalogs) == length(target_terms) )
PChidden_normal<- reference_PC_gen_matrix

for (j in 1:ncol(reference_PC_gen_matrix))
	PChidden_normal[,j] <- hidden_normalPC( ref_PC= reference_PC_gen_matrix[,j],
						target_metalog= target_metalogs[[j]],
						term= target_terms[[j]] )

PChidden_normal
}

#' @importFrom rmetalog qmetalog

pull2targetMetalog<- function(x_normal, target_metalog_distr, term) {
 u<- pnorm(x_normal, mean= 0, sd=1 )
qmetalog( target_metalog_distr, y= u, term=term) 
}

#' Get a set of cloned PCs transformed to target metalog distributions based on input normal PCs
#'
#' @details
#' Function will generate set of variables having same covariance as columns of \code{PChidden_normal}. Then, 
#' each such variable gets transformed to metalog distribution given by \code{target_metalogs[[i]]} 
#' (a \code{metalog} object).
#'
#' @param PC_hidden_normal A matrix of which covariance is computed and replicated in output data.
#' @param target_metalogs A list of \code{metalog} objects specifying target marginal distributions for generated variables. Must be of length equal to number of columns of \code{PC_hidden_normal}.
#' @param target_terms A vector of numbers of terms to use in applying each metalog transformation of the same length as \code{target_metalogs}.
#' @return A matrix of same dimension as \code{PC_hidden_normal} containing cloned version of input matrix which columns are uncorrelated with the input and are transformed to target metalog distributions.
#' @examples
#' data(brca)
#' data(brca_clusters)
#' lvl1<- initial_clusterNreconstruct(X= brca, X_variances=matrixStats::colVars(brca),
#'				        clustering_vector=brca_clusters)
#' PC_generator_matrix( lvl1$cluster_blockPCA) -> PCmat
#' meta_list<- list()
#' for (j in 1:ncol(PCmat))  rmetalog::metalog(PCmat[,j], term_limit=5, step_len=.01) -> meta_list[[j]]
#' HN_PCmat<- hidden_normalPC_matrix(PCmat, meta_list, rep(5, length(meta_list) ) )
#' cloned_metaPC<- cloned_metalogPC_matrix(HN_PCmat, meta_list, rep(5, length(meta_list)) )
#'  #note: no correlations between cloned and reference PCs
#' heatmap( abs(cor(cbind(PCmat, cloned_metaPC))), Colv=NA, Rowv=NA)
#' @export


cloned_metalogPC_matrix<- function( PChidden_normal, target_metalogs, target_terms){

stopifnot( ncol(PChidden_normal) == length(target_metalogs) )
stopifnot( length(target_metalogs) == length(target_terms) )
H_cov= cov(PChidden_normal)
H_clones<-multivariate_normal_cholesky(sigma=H_cov, n_samples=nrow(PChidden_normal) ) 
H_clones<- scale(H_clones, center=TRUE, scale=TRUE)
	colnames(H_clones)<- colnames(PChidden_normal)
PCmeta<- H_clones
for (j in 1:ncol(PCmeta))
	PCmeta[,j] = pull2targetMetalog(x_normal=H_clones[,j],target_metalog_distr=target_metalogs[[j]],
					term=target_terms[[j]]
					)
PCmeta 
}

add_noise2<- function( Xg, Eg_variances) {


for (j in 1:length(Eg_variances) )
	Xg[,j] = Xg[,j] + rnorm( nrow(Xg), 0, sqrt(Eg_variances[[j]]) )
return(Xg)
}

#' Simulate data using output of HCR and a set of synthetic principal components.
#'
#' @details This function does majority of the HCS method. It uses HCD of the input data \code{X} supplied in form of
#' \code{hcluster_blockPCA} and \code{hclustering_list} along with precomputed set of synthetic PCs of clusters
#' and subclusters in HCD (synthetic PCs are in \code{cloned_PCmatrix}. First two components can be produced by
#' functions \code{subclusterNreconstruct} or \code{initial_clusterNreconstruct}.
#' \code{cloned_PCmatrix} should contain as columns simulations of PCs of each cluster and subcluster,
#' independent from original data (can be produced by \code{cloned_normalPC_matrix} or \code{cloned_metalogPC_matrix}). 
#'
#' @param X An input dataset where rows are samples and columns variables
#' @param hcluster_blockPCA  

get_cloned_reconstruction<- function(X,hcluster_blockPCA, hclustering_list, cloned_PCmatrix, noise_variances=NULL, uncenter=TRUE) {

g<-get_max_g(hclustering_list)
print(g)
Xg<- X
Xg[,]=0  #placeholder for the reconstruction of gth order

for (g_j in 1:g) {
	get_partition_at_g( hclustering_list, g_j) -> g_j_clustering
	X_g_j <- X
	X_g_j[,]=0
	hcluster_blockPCA[ names(g_j_clustering) ]-> g_j_blocks
	names(g_j_blocks) <- names(g_j_clustering)
	for ( K in names(g_j_blocks) )
		g_j_blocks[[K]]$reconstruction <- cloned_PCmatrix[, colnames(cloned_PCmatrix)==K,drop=FALSE ] %*% g_j_blocks[[K]]$Vt
	allocate_blocks( where=X_g_j, blockwise_PCA= g_j_blocks )-> X_g_j
	Xg = Xg + X_g_j
}
if (!is.null(noise_variances)) Xg<- add_noise2(Xg, noise_variances)
if (uncenter) { colMeans(X)-> meansX
		for(j in 1:ncol(X)) Xg[,j] = Xg[,j] + meansX[[j]]
	      }
return(Xg)
}
