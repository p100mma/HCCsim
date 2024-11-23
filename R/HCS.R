

#' Extract matrix of principal components of all clusters of HCD
#' 
#' @param \code{hcluster_blockPCA} A list of blockwise PCA reconstruction of clusters and subclusters as obtained and \code{subclusterNreconstruct} function.
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
#' @seealso [initial_clusterNreconstruct()], [subclusterNreconstruct()] for how input data should look like
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
#' @seealso [MSE()]
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
#' @param n_samples Number of samples to generate, defaults to the number of rows of \code{reference_PC_gen_matrix}.
#' @return A matrix of with same number of columns as input one, containing synthetic PCs from normal distribution with identical covariance as reference ones but independent of original ones, with number of rows specified by \code{n_samples}.
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
#' @seealso [PC_generator_matrix()], [subclusterNreconstruct()], [initial_clusterNreconstruct()] for how input data should look like, [multivariate_normal_cholesky()] for data gen. method, [cloned_metalogPC_matrix()] for more elaborate simulation of PCs
#' @export


cloned_normalPC_matrix<- function(reference_PC_gen_matrix, n_samples= nrow(reference_PC_gen_matrix)) {
	PCsig<- cov(reference_PC_gen_matrix)
	PCnorm<-multivariate_normal_cholesky(sigma=PCsig, n_samples=n_samples ) 
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
#' First, each \code{reference_PC_gen_matrix[,i]} gets transformed by applying CDF of metalog [1] distribution
#' stored in \code{target_metalogs[i]} ( \code{metalog} object) using \code{target_terms[i]} terms.
#' Then, variable number \code{i} is pulled under standard normal distribution by applying inverse of its 
#' cumulative distribution function.
#'
#' @param reference_PC_gen_matrix Matrix where columns are PCs of clusters which will be transformed to normal distributions. 
#' @param target_metalogs A list of objects of class \code{metalog} of length equal to number of columns of \code{reference_PC_gen_matrix}.
#' @param target_terms A vector of numbers of terms to use in applying each metalog transformation of the same length as \code{target_metalogs}.
#' @return A matrix of the same dimensions as \code{reference_PC_gen_matrix}, containing transformed variables from input matrix, with same column names.
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
#' @seealso [PC_generator_matrix()] ,[subclusterNreconstruct()] ,[initial_clusterNreconstruct()] for how to generate input data for this function, [rmetalog::metalog()] for metalog distribution
#' @references [1] Keelin T.  The metalog distributions. Dec Anal  2016;13:243–77.
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
#' each such variable gets transformed to metalog distribution [1] given by \code{target_metalogs[[i]]} 
#' (a \code{metalog} object).
#'
#' @param PC_hidden_normal A matrix of which covariance is computed and replicated in output data. Column names of that matrix get copied to the output one.
#' @param target_metalogs A list of \code{metalog} objects specifying target marginal distributions for generated variables. Must be of length equal to number of columns of \code{PC_hidden_normal}.
#' @param target_terms A vector of numbers of terms to use in applying each metalog transformation of the same length as \code{target_metalogs}.
#' @param n_samples Number of samples to generate, defaults to the number of rows of \code{reference_PC_gen_matrix}.
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
#' @seealso [hidden_normalPC_matrix()], [PC_generator_matrix()] ,[subclusterNreconstruct()] ,[initial_clusterNreconstruct()] for how to generate input data for this function, [rmetalog::metalog()] for metalog distribution
#' @references [1] Keelin T.  The metalog distributions. Dec Anal  2016;13:243–77.
#' @export


cloned_metalogPC_matrix<- function( PChidden_normal, target_metalogs, target_terms, 
				    n_samples=nrow(PChidden_normal)){

stopifnot( ncol(PChidden_normal) == length(target_metalogs) )
stopifnot( length(target_metalogs) == length(target_terms) )
H_cov= cov(PChidden_normal)
H_clones<-multivariate_normal_cholesky(sigma=H_cov, n_samples=n_samples ) 
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
#' @param X An input dataset matrix where rows are samples and columns variables
#' @param hcluster_blockPCA Reconstructions of blocks defined by clusters from \code{hclustering_list}. Produced by \code{sublusterNreconstruct} or by \code{initial_clusterNreconstruct} or \code{blockwise_PCA_reconstruction}.
#' @param hclustering_list  A list encoding cluster membership, such as one produced by \code{subclusterNreconstruct} or \code{subdivide_cluster} or an instance of \code{HCCSim_clustering_list}.
#' @param cloned_PCmatrix Matrix of PCs simulating PCs in \code{hcluster_blockPCA}. Number of rows of that matrix determines number of samples to generate in output dataset. Can be generated by \code{cloned_normalPC_matrix} or \code{cloned_metalogPC_matrix} functions.
#' @param noise_variances A vector of variances of random noise to add to each synthetic variable. Has to be of length equal to number of columns of \code{X}. Such addition of random noise is performed if this argument is not NULL. Without adding random noise, generated data has lower variance and amplified correlations. 
#' @param uncenter If \code{TRUE}, then means of \code{X} will be added back to synthetic data. Otherwise, generated data is zero centered.
#' @return A matrix of dimension \code{nrow(cloned_PCmatrix)} x \code{ncol(X)}, containing synthetic data, generated by using PCs from \code{cloned_PCmatrix} and original coefficients of linear combinations of PCs from \code{hcluster_blockPCA}.
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
#'
#' ## 2 level HCS(n): normal distribution of PCs
#' PC_generator_matrix( lvl2$hcluster_blockPCA) -> PCmat2
#' cloned_PC_N<- cloned_normalPC_matrix(PCmat2)
#' #compute missing variance of 2lvl rec.
#' X_r2<- get_full_reconstruction(brca, lvl2$hcluster_blockPCA, lvl2$hclustering_list,
#' 				  add_noise=FALSE) #crucial arg to compute missing variance properly!
#' noise_variances2<- matrixStats::colVars(brca-X_r2)
#' X_s2<- get_cloned_dataset(brca, lvl2$hcluster_blockPCA, 
#'				lvl2$hclustering_list,
#'				cloned_PC_N,
#' 				noise_variances= noise_variances2)
#'
#' ## 1 lvl HCS(f): metalog distribution fitting on PCs 
#' PC_generator_matrix( lvl1$cluster_blockPCA) -> PCmat1
#' meta_list<- list()
#' for (j in 1:ncol(PCmat1))  rmetalog::metalog(PCmat1[,j], term_limit=5, step_len=.01) -> meta_list[[j]]
#' HN_PCmat<- hidden_normalPC_matrix(PCmat1, meta_list, rep(5, length(meta_list) ) )
#' ### note: we can generate any number of samples, say 1000
#' cloned_metaPC<- cloned_metalogPC_matrix(HN_PCmat, meta_list, rep(5, length(meta_list)),  n_samples= 1000 )
#' #compute missing variancce of 1lvl rec.
#' X_r1<- get_full_reconstruction(brca, lvl1$cluster_blockPCA, lvl1$clustering_list,
#' 				  add_noise=FALSE)
#' noise_variances1<- matrixStats::colVars(brca-X_r1)
#' X_s1<- get_cloned_dataset(brca, lvl1$cluster_blockPCA, 
#'				lvl1$clustering_list,
#'				cloned_metaPC,
#' 				noise_variances= noise_variances1)
#' @seealso [cloned_normalPC_matrix()], [cloned_metalogPC_matrix()] for how to generate synthetic PCs, [subclusterNreconstruct()] ,[initial_clusterNreconstruct()] for how to generate input data for this function, [get_full_reconstruction()] for reconstructing original data and producing output data correlated with original
#' @export



get_cloned_dataset<- function(X,hcluster_blockPCA, hclustering_list, cloned_PCmatrix, noise_variances=NULL, uncenter=TRUE) {

g<-get_max_lvl(hclustering_list)

Xg<- matrix(nrow=nrow(cloned_PCmatrix), ncol=ncol(X))
Xg[,]=0  #placeholder for the reconstruction of gth order

for (g_j in 1:g) {
	get_partition_at_g( hclustering_list, g_j) -> g_j_clustering
	X_g_j <- Xg
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

#' Simulate one (sub)module using input HCD of dataset and synthetic PCs
#'
#' @details This is modular version of [get_cloned_dataset()], see that function to get meaning of input data. 
#' 
#' Additional parameters `cluster_tuple` and `up_to` are passed to [extract_subclustering()]. 
#' They control which (sub)cluster to simulate and how many sublayers of it.
#'
#' If `add_noise=TRUE`, then missing variance in the module is complemented by random noise in the following manner:
#' The variance of each variable in the module will be eqial to
#' that variable in X, but variance contained in the module that carries
#'  contribution from PCs will only use PCs from cluster `cluster_tuple`
#' and `up_to` sublayers nontained in it.
#'
#' @param X An input dataset matrix where rows are samples and columns variables
#' @param cluster_tuple Name of the module (from `hclustering_list`) to simulate.
#' @param hcluster_blockPCA Reconstructions of blocks defined by clusters from \code{hclustering_list}. Produced by \code{sublusterNreconstruct} or by \code{initial_clusterNreconstruct} or \code{blockwise_PCA_reconstruction}.
#' @param hclustering_list  A list encoding cluster membership, such as one produced by \code{subclusterNreconstruct} or \code{subdivide_cluster} or an instance of \code{HCCSim_clustering_list}.
#' @param cloned_PCmatrix Matrix of PCs simulating PCs in \code{hcluster_blockPCA}. Number of rows of that matrix determines number of samples to generate in output dataset. Can be generated by \code{cloned_normalPC_matrix} or \code{cloned_metalogPC_matrix} functions.
#' @param noise_variances Variances of random noise to add to each variable in the module, must be of length `length(hclustering_list[[cluster_tuple]])`.
#' @param uncenter If \code{TRUE}, then means of \code{X} will be added back to synthetic data. Otherwise, generated data is zero centered.
#' @param up_to Number of sublayers of module `cluster_tuple` to include in the simulation. If `NULL`, are layers are included, if it is 0, only the part of the signal of the module ID-ed by `cluster_tuple` is included, if it exceeds the number of available sublayers of `cluster_tuple`, error is raised.
#' @return A matrix of dimension \code{nrow(cloned_PCmatrix)} x \code{length(hclustering_list[[cluster_tuple]])}, containing synthetic data according to `up_to` specification, generated by using PCs from \code{cloned_PCmatrix} and original coefficients of linear combinations of PCs from \code{hcluster_blockPCA}.
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
#'
#' ## 2 level HCS(n): normal distribution of PCs
#' PC_generator_matrix( lvl2$hcluster_blockPCA) -> PCmat2
#' cloned_PC_N<- cloned_normalPC_matrix(PCmat2)
#' #compute missing variance of 2lvl rec.
#' X_r2<- get_full_reconstruction(brca, lvl2$hcluster_blockPCA, lvl2$hclustering_list,
#' 				  add_noise=FALSE) #crucial arg to compute missing variance properly!
#' noise_variances2<- matrixStats::colVars(brca-X_r2)
#' #First - full simulation. Then - do modular version and compare results for the module, should be similar
#' X_s2<- get_cloned_dataset(brca, lvl2$hcluster_blockPCA, 
#'				lvl2$hclustering_list,
#'				cloned_PC_N,
#' 				noise_variances= noise_variances2)
#' module_2<-simulate_module(X=brca, cluster_tuple='2',
#' 			     hcluster_blockPCA= lvl2$hcluster_blockPCA,
#'			     hclustering_list=lvl2$hclustering_list,
#'			     cloned_PCmatrix=cloned_PC_N,
#'			     noise_variances= noise_variances2[ lvl2$hclustering_list[['2']] ])
#' c1<- corfast(X_s2[, lvl2$hclustering_list[['2']]])	
#' c2<- corfast(module_2)
#' MSE(c1,c2)
#' cc1<- ccfast(c1^2)
#' cc2<- ccfast(c2^2)
#' par(mfrow=c(1,3))
#' hist(cc1, main="local cluster coef., in module 2 from full sim")
#' plot(cc1,cc2, main="cc module 2, full sim vs modular")	
#' hist(cc2, main="local cluster coef., in module 2 from full modular ver.")
#' @export
simulate_module<- function(X, cluster_tuple,
 			   hcluster_blockPCA,
			   hclustering_list,
			   cloned_PCmatrix,
			   noise_variances=NULL,
			   uncenter=TRUE,
			   up_to=NULL
			   ) {
	cloned_PCmatrix<- cloned_PCmatrix[sample(1:nrow(cloned_PCmatrix), nrow(cloned_PCmatrix), replace=FALSE),]
	X_s= matrix(nrow=nrow(cloned_PCmatrix), ncol= ncol(X) ) 
	X_s[,]=0
	colnames(X_s)<- colnames(X)
	subclusters<-extract_subclustering( cluster_tuple, hclustering_list,up_to)	
	subclusters_PCA<- hcluster_blockPCA[ names(subclusters) ]
	get_max_lvl(subclusters) -> end_lvl
	# HCR.R for p.() definition
	length(p.(cluster_tuple)) -> start_lvl
	for (layer_number in start_lvl:end_lvl) {
	   # since cluster placements are in temrs of X,
	   # we will put the subcluster here and transfert it to mod then
	   X_aux= matrix(nrow=nrow(cloned_PCmatrix), ncol= ncol(X) ) 
	   X_aux[,]=0
	   subclusters_thisLayer<- get_partition_at_g(subclusters, layer_number)
	   print(names(subclusters_thisLayer))
	   subclusters_PCA_thisLayer<- subclusters_PCA[ names(subclusters_thisLayer) ]
	   print(names(subclusters_PCA_thisLayer))
	   for (K in names(subclusters_PCA_thisLayer) )
	  	subclusters_PCA_thisLayer[[K]]$reconstruction <- cloned_PCmatrix[, colnames(cloned_PCmatrix)==K,drop=FALSE] %*% subclusters_PCA_thisLayer[[K]]$Vt
	    allocate_blocks(where=X_aux, blockwise_PCA=subclusters_PCA_thisLayer) -> X_aux
	    X_s = X_s + X_aux 
	   }
	# now we take just the subcluster
	mod=X_s[, subclusters[[cluster_tuple]], drop=FALSE ] 
	if (uncenter | (!is.null(noise_variances)))
		X_mod<- X[, subclusters[[cluster_tuple]], drop=FALSE ] 
	if (!is.null(noise_variances)) mod<-add_noise2(mod,noise_variances)
	if (uncenter) { colMeans(X_mod)-> meansX
		for(j in 1:ncol(X_mod)) mod[,j] = mod[,j] + meansX[[j]]
	      }
	return(mod)
}

