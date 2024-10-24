

#' Perform k-th order PCA reconstruction of correlated data with clustered variables
#'
#' In a nutshell, function performs PCA separately in each group of columns of \code{A} defined by \code{clustering_list} using at most \code{k} PCs per group of columns, trying to explain a fraction \code{f} of total variance inside each group and uses the result to reconstruct \code{A} using those clusterwise PCs.
#'
#' @param A Data matrix where columns are variables and rows are samples.
#' @param A_variances Precomputed variances of columns of \code{A}. If not computed already, set this to \code{colVars(A)}
#' @param clustering_list A named list where names are cluster labels (integers converted to strings) and elemnts are vectors of indexes of nodes (columns of \code{A}) in each cluster, as produced by \code{clustering_vector2list} or outlined in \code{HCCSim_clustering_list}.
#' @param f Target fraction of variance to explain, a real number in \code{(0,1]}.
#' @param k Maximum number of PCs to use per cluster in \code{clustering_list}.
#' @return A \code{list} with one element per cluster in \code{clustering_list}. Each such element contains the reconstruction of columns of \code{A} in that cluster and additional information on the result, i.e. each element is a list with the following elements: 
#' \itemize{
#' \item \code{location} A vector of indexes of reconstructed columns in the block as placed in \code{A}.
#' \item \code{total_var} Real number, total variance of the block.
#' \item \code{pc_var} Variances of all PCs of each block.
#' \item \code{cumul_var} Vector of cumulative variances of all PCs.
#' \item \code{k_used} Number of PCs used (can be less than \code{k} )
#' \item \code{k_f} Number of PCs one would have to use to explain target fraction \code{f} of total variance.
#' \item \code{f} Same as input argument.
#' \item \code{saturated} True if total variance was fully explained with \code{k_used} PCs.
#' \item \code{var_explained} How much variance in total was explained.
#' \item \code{PC} Matrix of PCs used for reoconstruction, has \code{k2use} columns, one per each PC.
#' \item \code{Vt} Coefficients of PCs used for reconstruction.
#' \item \code{reconstruction} Reconstructed data using \code{PC} and \code{Vt}.
#' 	   }
#' @importFrom matrixStats colVars
#' @examples
#' data(brca)
#' data(brca_clusters)
#' cl_list<- clustering_vector2list( brca_clusters ) 
#' blocks_brca<- blockwise_PCA_reconstruction(A=brca, A_variances= matrixStats::colVars(brca), clustering_list= cl_list, f=0.6, k=3)
#' @export

blockwise_PCA_reconstruction<- function(A,A_variances, clustering_list, f, k){
   lapply(clustering_list, function(C_i)
        {
            A.cl<- A[, C_i, drop=FALSE]
            A.cl<-scale(A.cl, center=TRUE, scale=FALSE)
	    sum(A_variances[ C_i])-> var.total
	    svd(A.cl)-> svd_A.cl
	    PC_A.cl<- with(svd_A.cl,  u %*% diag(d, ncol=length(d), nrow=length(d)) )
	    pc.vars<-colVars(PC_A.cl)
            cumsum(pc.vars)-> pc.vars.cumul
	    if( length(C_i)==1 ) {k_i =1; k2use=1} else {
            k_i <-min(which(pc.vars.cumul/var.total >= f  ))
	    k2use<- min(k_i, k) }
	    lc_coefs_ALL<- t(svd_A.cl$v)
	    lc_coefs<- lc_coefs_ALL[1:k2use,, drop=FALSE]
	    PC_rec<- PC_A.cl[, 1: k2use, drop=FALSE ]
	    PC_rec %*%  lc_coefs -> A.cl_rc
            list(location=C_i,
                 total_var=var.total,
		 pc_var=pc.vars,
		 cumul_var=pc.vars.cumul,
                 k_used=k2use,
		 k_f=k_i,
		 f=f,
		 saturated= k_i <= k2use,
                 var_explained= pc.vars.cumul[[k2use]],
                 PC= PC_rec,
                 Vt=lc_coefs,
	#	 all_Vt=lc_coefs_ALL,
		 reconstruction=A.cl_rc
                 )
	})-> blocks_list; names(blocks_list)<- names(clustering_list)
		return(blocks_list)
								}

#' Place blockwise PCA reconstructions in an input container matrix
#'
#' @param where Input matrix where rows are samples and columns variables. It should have the same dimensions as matrix used as input to \code{blockwise_PCA_reconstruction} function used to produce \code{blockwise_PCA} argument.
#' @param blockwise_PCA Result of \code{blockwise_PCA_reconstruction} function ot use.
#' @return An updated matrix \code{where} with columns corresponding to each cluster were replaced by their recondstructions from \code{blockwise_PCA} objects.
#' @examples
#' data(brca)
#' data(brca_clusters)
#' cl_list<- clustering_vector2list( brca_clusters ) 
#' blocks_brca<- blockwise_PCA_reconstruction(A=brca, A_variances= matrixStats::colVars(brca), clustering_list= cl_list, f=0.6, k=3)
#' brca_rec<- brca
#' brca_rec[,]=0
#' brca_rec= allocate_blocks( brca_rec, blocks_brca)
#' @export
allocate_blocks<-function( where, blockwise_PCA){

    for (block in blockwise_PCA){
        where[, block$location] <- block$reconstruction
			    }
where
						}

#' Recalculate fraction of variance to explain taking the clusterless variables into the account
#'
#' @param clustering_vector An integer vector for which \code{clustering_vector[i]} gives a label of the group to which node \code{i} belongs. Label of \code{0} should encode clusterless nodes.
#' @param f_E Target fraction of total variance in the dataset to explain.
#' @param X_variances Precalculated variances of each variable in the dataset, should be the same length and ordering as in \code{clustering_vector}. If not calculated already, set this to \code{colVars(X)} where \code{X} is matrix where rows are samples and columns variables.
#' @return A number in \code{(0,1]}, recalulated \code{f_E} such that the recalculated fraction takes into the account only the variables in proper clusters (label \code{!=0} ).
#' @examples
#' data(brca)
#' data(brca_clusters)
#' calculate_f( zeroOutSmall(brca_clusters, 30), f_E=0.7, X_variances= matrixStats::colVars(brca) )
#' @export

calculate_f<- function( clustering_vector, f_E, X_variances){
 v_total<- sum(X_variances)
 v_noiseless<- sum(X_variances[ clustering_vector!= 0 ])
 v_target<- f_E * v_total
 if (v_noiseless < v_target) {
	warning("initial target fraction f_E cannot be reached by using any fraction of noiseless clustering, returning NA")
	return(NA)
    }
 stopifnot(v_noiseless>0)
 return( v_target/ v_noiseless)	
}


#######EASY PART- one level reconstruction ##########################


#' Divide variables dataset based on input clustering and reconstruct clusters using at most k principal components
#' 
#' @param X A matrix containing dataset where columns are variables and rows are samples.
#' @param X_variances Precomputed variances of \code{X}.
#' @param clustering_vector A vector encoding cluster membership according to convention outlined in \code{HCCsim_clustering_vector} class. Either this or \code{clustering_list} argument shall be non null.
#' @param clustering_list  A list encoding cluster membership according to convention outlined in \code{HCCsim_clustering_list} class. Either this or \code{clustering_list} argument shall be non null.
#' @param min_cl_size If given, nodes in clusters of size less than \code{min_cl_size} get classified as clusterless nodes.
#' @param f_E A target fraction of variance to explain. Gets recalculated to match the corresponding fraction of total variance of proper clusters (without clusterless variables).
#' @param k Maximum number of principal components to use for reconstructing each cluster.
#' @return A \code{list} with following elements:
#' \itemize{
#' \item \code{clustering_list} A \code{clustering_list} representing input clustering.
#' \item \code{f} Recalculated \code{f_E} giving the fraction of variance of proper clusters to explain (excluding "noise" clusters (\code{C_0}))
#' \item \code{cluster_blockPCA} A \code{k}-th order clustering based reconstruction of \code{X}, based on \code{clustering_list} component. It is output of \code{blockwise_PCA_reconstruction} funtion.
#' }
#' @examples
#' data(brca)
#' data(brca_clusters)
#' result<- initial_clusterNreconstruct(X= brca, X_variances=matrixStats::colVars(brca),
#'				        clustering_vector=brca_clusters)
#' lapply(result$clustering_list, length)
#' length(result$cluster_blockPCA)
#' @export

initial_clusterNreconstruct<- function(X, 
				       X_variances,
				       clustering_vector=NULL,
				       clustering_list=NULL,
				       min_cl_size=min(30,ncol(X)),
				       f_E=0.6,
				       k=5) 
	{
	##correctness check for input 
	stopifnot( length(f_E)==1)
	stopifnot( (f_E > 0 ) && (f_E <= 1) )
	stopifnot(ncol(X)==length(X_variances))	
	if (is.null(clustering_vector) && is.null(clustering_list))
		stop("One of: clustering_vector, clustering_list must be not null")
	##based on given input, create another one needed
	if (is.null(clustering_list))
		{
		 stopifnot(ncol(X)==length(clustering_vector))
		 clustering_list<- clustering_vector2list(clustering_vector)
		}
	if (is.null(clustering_vector))
		{
		stopifnot( all(unlist(clustering_list) %in% 1:ncol(X) ) )
		 clustering_vector<- clustering_list2vector(clustering_list)
		}	
	##separate C_0 if min_cl_size given
	if (!is.null(min_cl_size))
		{
		stopifnot(min_cl_size <= ncol(X) )
		clustering_vector<-zeroOutSmall(clustering_vector, min_cl_size)
		#relabel proper clusters to 1:N_C and update clustering_list
		tidyUpLabels(clustering_vector)-> clustering_vector
		 clustering_list<- clustering_vector2list(clustering_vector)
		}
	f=calculate_f( clustering_vector, f_E, X_variances)

	# if calculation of 'f' failed:	
	if (is.na(f)) {warning(" stopped because provided clustering leaves too much variance in the noise part");
	       return( list(clustering_list=clustering_list,
			    f=f,
			    cluster_blockPCA=NA) ) }
		
        # result if everything went OK	
     list(clustering_list=clustering_list, f=f,
     cluster_blockPCA=blockwise_PCA_reconstruction(X,
						  X_variances, 
						  clustering_list,
						  f=f, k=k)
	)
} 
	
	
	
#######GENERALIZATION TO HIERARCHICAL STRUCTURE #########################

#utility function, indexes 1,2,3 to tuple '1,2,3'
#idempotent, i.e .p(.p(1,2,3)) = .p('1,2,3') = '1,2,3'
#clusters will be indexed by such strings of tuples
.p<- function(...) paste0(c(...),collapse=',')

#divide hclustering_list or
#initialize hclustering_list by subdividing clustering_list



#' Add subclusters to list describing hierarchical clustering with correct indexing according to tuple convention
#'
#' @param hclustering_list A list of clusters and sub-clusters indexed by strings of integer tuples (like '1','2','1,1','1,2','2,1','2,2' etc.). Each cluster is a vector of integers dentoing nodes belonging to that cluster. Can be an object of class \code{HCCsim_clustering_list}.
#' @param cluster_tuple A vector of indexes describing the indexing tuple or a string (c(1,2) or '1,2')
#' @param subdivision A grouping of indexes inside cluster indexed by \code{cluster_tuple} into subgroups in form of \code{clustering_list}.
#' @return An updated \code{hclustering_list} with additional elements denoting subclusters of cluster indexed by \code{cluster_tuple}, indexed by the tuple convention ( e.g. if \code{cluster_tuple='1,2'} then added subclusters from subdivision with clusters '1','2','3' are '1,2,1', '1,2,2', '1,2,3'). Integer suffixes are inferred from the order of the clusters in \code{subdivision}.
#' @examples
#' base_cl<- HCCSim_clustering_list(list( `1`= c(5,6,7,8), `2`= c(9,4,3), `3`=c(1,2,10,11) ),
#' 				    domain_size= 15,   #4 clusterless nodes implied
#' 				    params=list(base_param=1) )
#' subdivision1= HCCSim_clustering_list(list(c(6,8), c(7,5) ), params=list(sub1_param=2))
#' subdivision2= HCCSim_clustering_list(list(c(1,10), c(2,11) ), params=list(sub2_param=3))
#' subdivide_cluster( base_cl, "1", subdivision1)-> hcl
#' print(hcl)
#' subdivide_cluster( hcl, "3", subdivision2)-> hcl
#' print(hcl)
#' @export



subdivide_cluster<- function( hclustering_list,
			      cluster_tuple,
                                subdivision){

if (length(hclustering_list)==0) stop(" hclustering_list must be nonempty") 
#now cluster_tuple is surely a string vector with 1 element
cluster_tuple<- .p( cluster_tuple)
if (!(cluster_tuple %in% names(hclustering_list))) 
	stop("cluster to divide is not an element of hclustering_list!")
 stopifnot( all ( unlist(subdivision) %in% hclustering_list[[ cluster_tuple ]] ) )
 stopifnot( all ( unlist(lapply(subdivision, length)) > 0 ) )
 stopifnot(  length( unlist( unique(subdivision)) ) == length(unlist(subdivision)))
 for (i in seq_along(subdivision))
	hclustering_list [[ .p(cluster_tuple, i) ]]<- subdivision[[i]]
attr(subdivision,'params')-> s_params
names(s_params)<- paste0( cluster_tuple,'_',  names(s_params))
attr(hclustering_list, 'params') <- c(attr(hclustering_list, 'params'), s_params)
 hclustering_list
}

#inverse of .p
p.<-function(str_tuple) unlist(strsplit(str_tuple, split=','))


#' Extract number of levels of hierarchy from \code{hclustering_list}

#' @param hclustering_list A list of clusters and sub-clusters indexed by strings of integer tuples (like '1','2','1,1','1,2','2,1','2,2' etc.). Each cluster is a vector of integers dentoing nodes belonging to that cluster. Designed to work with output of \code{subdivide_cluster} function or instance of \code{HCCsim_clustering_list} class.
#' @return An integer giving the number of levels of hierarchy in \code{hclustering_list}
#' @examples
#' base_cl<- HCCSim_clustering_list(list( `1`= c(5,6,7,8), `2`= c(9,4,3), `3`=c(1,2,10,11) ),
#' 				    domain_size= 15,   #4 clusterless nodes implied
#' 				    params=list(base_param=1) )
#' subdivision1= HCCSim_clustering_list(list(c(6,8), c(7,5) ), params=list(sub1_param=2))
#' subdivision2= HCCSim_clustering_list(list(c(1,10), c(2,11) ), params=list(sub2_param=3))
#' subdivide_cluster( base_cl, '1', subdivision1)-> hcl
#' subdivide_cluster( hcl, '3', subdivision2)-> hcl
#' get_max_lvl(hcl)
#' @export

get_max_lvl<- function(hclustering_list){
if (length(hclustering_list)==0) stop(" hclustering_list must be nonempty") 
names(hclustering_list)-> all_tuples
 tuple_lengths<-unlist( lapply(all_tuples, function(tpl) length(p.(tpl)) ))
 max(tuple_lengths)
}

#' Extract \code{clustering_list} describing division into the clusters at specified lvl of hierarchy of \code{hclustering_list}
#' 
#' @param hclustering_list A list of clusters and sub-clusters indexed by strings of integer tuples (like '1','2','1,1','1,2','2,1','2,2' etc.). Each cluster is a vector of integers denoting nodes belonging to that cluster. Designed to work for output of \code{subdivide_cluster} function or an instance of \code{HCCsim_clustering_list} class.
#' @param g An integer specifying lvl of hierarchy (number of subdivision step) to extract from \code{hclustering_list} from which clusters to extract.
#' @return A \code{clustering_list} describing division into the clusters at lvl \code{g} of hierarchy of \code{hclustering_list} 
#' @examples
#' base_cl<- HCCSim_clustering_list(list( `1`= c(5,6,7,8), `2`= c(9,4,3), `3`=c(1,2,10,11) ),
#' 				    domain_size= 15,   #4 clusterless nodes implied
#' 				    params=list(base_param=1) )
#' subdivision1= HCCSim_clustering_list(list(c(6,8), c(7,5) ), params=list(sub1_param=2))
#' subdivision2= HCCSim_clustering_list(list(c(1,10), c(2,11) ), params=list(sub2_param=3))
#' subdivide_cluster( base_cl, '1', subdivision1)-> hcl
#' subdivide_cluster( hcl, '3', subdivision2)-> hcl
#' get_partition_at_g(hcl,1)
#' get_partition_at_g(hcl,2)
#' @export

get_partition_at_g<- function( hclustering_list,
			      g){
 stopifnot(g>0)
if (length(hclustering_list)==0) stop(" hclustering_list must be nonempty") 
names(hclustering_list)-> all_tuples
 tuple_lengths<-unlist( lapply(all_tuples, function(tpl) length(p.(tpl)) ))
 stopifnot( g <= max(tuple_lengths))
 g_tuples<- all_tuples[ tuple_lengths == g ] 
 g_P<-hclustering_list [ g_tuples ] 
 attr(g_P,'params')<- attr(hclustering_list, 'params')
 g_P
}

#' Calculate entropy of a division of a set into groups
#' 
#' @param clustering_list A list where elements are atomic vectors describing separate clusters.  Number of elements in each vector is the number of objects inside that cluster.
#' @return A real number giving the entropy of proportions of groups in the input.
#' @examples
#' base_cl<- HCCSim_clustering_list(list( `1`= c(5,6,7,8), `2`= c(9,4,3), `3`=c(1,2,10,11) ),
#' 				    params=list(base_param=1) )
#' normalized_cl_entropy(base_cl)
#' @export

normalized_cl_entropy<- function( clustering_list ) {
 unlist(lapply(clustering_list, length))-> n_c
 p_c<- n_c/sum(n_c)
 length(n_c)-> n_groups
 -sum( p_c * log( p_c ) )/(log(n_groups))
}

#' Normalized version of hierarchical clustering suited for use with HCR
#'
#' Clustering functions used with HCR must take as an input similarity matrix and return vector of cluster labels. This function adapts \code{hclust} from base R to this convention.
#'
#' @param similarity_matrix A symmetric matrix of node similarities.
#' @param n_group A number of groups in the output clustering.
#' @param max_sim Maximum value of similarities in \code{similarity_matrix}.
#' @param method A method for cluster to cluster distance calculation in \code{hclust}.
#' @return An integer vector of cluster labels of nodes.
#' @examples
#' data(brca)
#' table(similarity_based_hclust( cor(brca)^2, 7))
#' @export 

similarity_based_hclust<- function( similarity_matrix, n_group, max_sim=1, method="complete") {
hclust( as.dist( max_sim - similarity_matrix), method=method)-> hcl_object
clvec<-cutree(hcl_object, k=n_group)
attr(clvec, 'params')=list(n_group=n_group)
attr(clvec, 'domain_size')= ncol(similarity_matrix)
clvec
}

#' Subdivide clusters in the innermost layer of hierarchical clustering and reconstruct the residuals in the subdivided clusters
#' 
#' @param X A matrix containing dataset where columns are variables and rows are samples.
#' @param X_variances Precomputed variances of \code{X}.
#' @param hclustering_list  A list encoding cluster membership, such as one produced by this function or \code{subdivide_cluster} or an instance of \code{HCCSim_clustering_list}. Clusters at the lowest level will be subdivided and reconstructed.
#' @param hcluster_PCA_blocks Corresponding reconstructions of blocks defined by clusters from \code{hclustering_list}. Produced by this function or by \code{initial_clusterNreconstruct} or \code{blockwise_PCA_reconstruction}.
#' @param k Maximum number of principal components to use for reconstructing each cluster.
#' @param f Target fraction of variance to explain (of the proper clusters).
#' @param clfun2 Clustering function to group residuals of each innermost cluster of \code{hclustering_list} based on magnitude of their correlations. Must receive similarity matrix as a first argument. Best clustering from several different settings will be picked by maximizing \code{normalized_cl_entropy()} value.
#' @param clfun2OtherArgs_constant List of named arguments of \code{clfun2} apart from similarity matrix, which shall be kept constant across different candidate settings to test.
#' @param clfun2OtherArgs_ranges A list of named vectors of parameters to test. Each named vector in the list shall have the same length. Names of the vectors should correspond to values of named arguments of \code{clfun2}. At step \code{i}, settings of \code{clfun2} are tested defined by \code{i}-th entry of each named vector in that list. Finally, for each cluster, setting maximizing value of \code{normalized_cl_entropy()} function will be chosen as definitive subdivision to use.
#' @param S_power A value of exponent which to use for transforming correlation of residuals to similarity by \code{abs(cor_matrix)^S_power}.
#' @param corfun A function to compute correlation of residuals of each cluster. Must accept matrix as its only argument (where columns are variables and rows are samples.
#' @return A \code{list} with following elements:
#' \itemize{
#' \item \code{hclustering_list} An expanded list representing resulting hierarchical clustering, it contains additional clusters that resulted from subdivision.
#' \item \code{hcluster_blockPCA} A \code{k}-th order clustering based reconstruction of residuals of reconstructions of previous level (determined by input \code{hcluster_blockPCA}. It is expanded by including reconstructions of subclusters added to expanded \code{hclustering_list}.
#'  }
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
#' @export

subclusterNreconstruct<- function(X,
				X_variances, hclustering_list, hcluster_PCA_blocks, k=5, f=0.6,
				clfun2,
				 clfun2OtherArgs_constant=NULL,
				  clfun2OtherArgs_ranges=NULL,
				  S_power=2,
				  corfun=cor){
Xg<-X
Xg[,]=0 # will get filled up by reconstruction up to current 'g' (max lvl)
g<-get_max_lvl(hclustering_list)
#Filling up Xg with reconstructions of each lvl...
for (g_j in 1:g) {
	get_partition_at_g( hclustering_list, g_j) -> g_j_clustering
	X_g_j <- X
	X_g_j[,]=0
	allocate_blocks( where=X_g_j, blockwise_PCA= hcluster_PCA_blocks [ names(g_j_clustering) ] )-> X_g_j
	Xg = Xg + X_g_j
} #...done
Xg_variances<- colVars(Xg)
# residuals of g-th order:
Eg= X - Xg 
Eg_variances<- colVars(Eg)
g_partition<- get_partition_at_g( hclustering_list, g)
g_blockPCA<- hcluster_PCA_blocks[ names(g_partition) ] 
for (bl in 1:length(g_blockPCA)) #if 'exhausted' status was not set for a block, declare it to be FALSE
	if (is.null(g_blockPCA[[bl]]$exhausted)) g_blockPCA[[bl]]$exhausted = FALSE

#this will (possibly) be expanded by subdivisions of existing clusters in g_partition
hcl_expanded<- hclustering_list
hcl_expanded_blockPCA<- hcluster_PCA_blocks

#get clusters which have some variance left and are big enough

##these were too small
exhausted_clusters<- names(g_blockPCA) [which( unlist(lapply(g_blockPCA, function(x) x$exhausted) )) ]
##remove above from these
unsaturated_tuples<- names(g_blockPCA) [which( unlist(lapply(g_blockPCA, function(x) !x$saturated) )) ]
unsaturated_tuples<- setdiff ( unsaturated_tuples, exhausted_clusters)    
#... now these tuples can be used to index hclustering_list

if (!length(unsaturated_tuples)) stop(' all clusters are either exhausted, or unsaturated. stopping')

# loop through clusters
#   calculate how much of variance of each cluster is to be explained still
#    then subdivide residuals in each cluster into subclusters 
#    by optimizing entropy of the division
#    and reconstruct it using at most 'k' components
#    if the size of cluster to subdivide is less than 'k' mark it as exhausted
#    and input as last lvl  the whole cluster itself
for (K in unsaturated_tuples) {
	C_K<- g_partition[[K]]
        left2explain= f* sum(X_variances[ C_K ] ) - sum(Xg_variances[ C_K] )
	#if estimate reaches 0 skip cluster C_K and declare it saturated
	if (left2explain <=0)
		{
		 hcl_expanded_blockPCA[[ K ]]$saturated= TRUE
	} else  {
		#one last check: is C_K big enough?
		if (length(C_K) <= k) 
		{
		 hcl_expanded_blockPCA[[ K ]]$exhausted= TRUE
		} else {
		# variance left2explain > 0 and C_K big enough so 
		# calculate similarity graph
		# of residuals & split
		S_Eg.C_K<-  abs (corfun( Eg [, C_K ])	)^S_power
		# choose best clustering among those produced by 
		# different values of clfun2OtherArgs_ranges
		clf2Args<- c(list( S_Eg.C_K), clfun2OtherArgs_constant )
		settings_scores<-vector()
		candidate_cl_lists<-list()
		settings_m<- length(clfun2OtherArgs_ranges[[1]])
		for (i_settings in 1:settings_m) {
			for (pname in names(clfun2OtherArgs_ranges))
				clf2Args[[pname]] = clfun2OtherArgs_ranges[[pname]][[i_settings]]	
			clustering_vector2list(do.call(clfun2, clf2Args))-> candidate_cl
			can_attrs<- attributes(candidate_cl)
			lapply(candidate_cl, function(C_K.j)  C_K[ C_K.j ] ) -> candidate_cl #translate local to global indexes
			attributes(candidate_cl)<- can_attrs
			candidate_cl_lists[[ i_settings ]] = candidate_cl
			settings_scores [[ i_settings ]] =  normalized_cl_entropy( candidate_cl )     # "best" percollating clustering is the one which proportions are closest
			     #... to an even split
			}	
		 C_K_subdivision<- candidate_cl_lists [[ which.max( settings_scores) ]]	
	noiseless_C_K_var<- sum(colVars ( Eg[,unlist(C_K_subdivision), drop=FALSE]))
	# rename each cluster 'j' in subdivision from 'j' to 'i1,i2..ig,j where K= i1,i2...ig' 
	names(C_K_subdivision) <- unlist(lapply( names(C_K_subdivision), function(j) .p( K, j ) ) )
	subdivide_cluster( hclustering_list=hcl_expanded,
			      cluster_tuple=K,
                                subdivision= C_K_subdivision) -> hcl_expanded #add subclusters to the hclust list
	f_K = left2explain/noiseless_C_K_var #target fraction of variance for the reconstruction of Eg in C_K
	if (left2explain > noiseless_C_K_var) stop('something is wrong, left over variance of residuals is smaller than the left2explain value')
	blockwise_PCA_reconstruction(Eg,Eg_variances, clustering_list= C_K_subdivision,
				      f=f_K, 
					k=k) -> C_K.blockPCA	
 	hcl_expanded_blockPCA<- c( hcl_expanded_blockPCA, C_K.blockPCA)	
	}#else of size of C_K<=0 
	}#else of left2explain <=0 END 
	} #C_K cluster subdivision END
     return( list( hclustering_list= hcl_expanded, #touched_tuples=touched_tuples,
	      hcluster_blockPCA= hcl_expanded_blockPCA) )
}

#' Complement the missing variances of reconstructed variables by appropriate amount of independent white noise
#'
#' Add \code{E_i} of noise from normal distribution to each reconstructed variable \code{Xg_i} such that \code{var(Xg_i)==var(X_i)} where \code{X_i} is a reference variable.
#'
#' @param Xg l
#' @param X l
#' @return An updated matrix \code{Xg} with random noise added to each of the columns.
#' @examples
#' data(brca)
#' data(brca_clusters)
#' lvl1<- initial_clusterNreconstruct(X= brca, X_variances=matrixStats::colVars(brca),
#'				        clustering_vector=brca_clusters)
#' rec<- brca
#' rec[,]=0
#' rec= allocate_blocks( rec, lvl1$cluster_blockPCA)
#' max( matrixStats::colVars(brca) - matrixStats::colVars(rec) ) #missing variances
#' rec=add_noise(rec, brca)
#' max( matrixStats::colVars(brca) - matrixStats::colVars(rec) ) #after complementing with random noise
#' @export

add_noise<- function( Xg, X) {

Eg <- X - Xg
Eg_variances<- apply(Eg,2,var)

for (j in 1:length(Eg_variances) )
	Xg[,j] = Xg[,j] + rnorm( nrow(Eg), 0, sqrt(Eg_variances[[j]]) )
return(Xg)
}




#' Perform the full reconstruction from HCD of a dataset
#'
#' Given a HCD (hierarchical clustering decomposion of correlation structure) of a dataset, generate it's reconstruction (HCR) based on that HCD.
#'
#' @param X A matrix containing dataset where columns are variables and rows are samples.
#' @param hcluster_blockPCA Reconstructions of blocks defined by clusters from \code{hclustering_list}. Produced by \code{sublusterNreconstruct} or by \code{initial_clusterNreconstruct} or \code{blockwise_PCA_reconstruction}.
#' @param hclustering_list  A list encoding cluster membership, such as one produced by \code{subclusterNreconstruct} or \code{subdivide_cluster} or an instance of \code{HCCSim_clustering_list}.
#' @param add_noise If \code{TRUE}, missing variance of reconstructed variables will be complemented by an adequate amount of random noise.
#' @param uncenter If \code{TRUE}, the means of variables in \code{X} will be added to corresponding variables in reconstruction.
#' @return A matrix \code{Xg} which is the \code{g}-th order hierarchical clustering reconstruction of correlated variables in \code{X}, where \code{g} is the number of levels in \code{hclustering_list} and \code{hcluster_blockPCA}.
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
#' X_2<- get_full_reconstruction(brca, lvl2$hcluster_blockPCA, lvl2$hclustering_list)
#' @export

get_full_reconstruction<- function(X,hcluster_blockPCA, hclustering_list, add_noise=TRUE, uncenter=TRUE) {

g<-get_max_lvl(hclustering_list)
Xg<- X
Xg[,]=0  #placeholder for the reconstruction of gth order

for (g_j in 1:g) {
	get_partition_at_g( hclustering_list, g_j) -> g_j_clustering
	X_g_j <- X
	X_g_j[,]=0
	allocate_blocks( where=X_g_j, blockwise_PCA= hcluster_blockPCA [ names(g_j_clustering) ] )-> X_g_j
	Xg = Xg + X_g_j
}
if (add_noise) Xg<- add_noise(Xg, X)
if (uncenter) { colMeans(X)-> meansX
		for(j in 1:ncol(X)) Xg[,j] = Xg[,j] + meansX[[j]]
	      }
return(Xg)
}


