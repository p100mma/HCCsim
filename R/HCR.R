

#' Perform k-th order PCA reconstruction of correlated data with clustered variables
#'
#' In a nutshell, function performs PCA separately in each group of columns of \code{A} defined by \code{clustering_list} using at most \code{k} PCs per group of columns, trying to explain a fraction \code{f} of total variance inside each group and uses the result to reconstruct \code{A} using those clusterwise PCs.
#'
#' @param A Data matrix where columns are variables and rows are samples.
#' @param A_variances Precomputed variances of columns of \code{A}. If not computed already, set this to \code{colVars(A)}
:q
:q
#' @param clustering_list A named list where names are cluster labels (integers converted to strings) and elemnts are vectors of indexes of nodes (columns of \code{A}) in each cluster, as produced by \code{clustering_vector2list} or outlined in \code{HCCSim_clustering_list}.
#' @param f Target fraction of variance to explain, a real number in \code{(0,1]}.
#' @param k Maximum number of PCs to use per cluster in \code{clustering_list}.
#' @return A \code{list} with one element per cluster in \code{clustering_list}. Each such element contains the reconstruction of columns of \code{A} in that cluster and additional information on the result, i.e. each element is a list with the following elements: 
#' \itemize{
#' \item \code{location} A vector of indexes of reconstructed columns in the block as placed in \code{A}.
#' \item \code{total_var} Real number, total variance of the block.
#' \item \code{pc_var}
#' \item \code{cumul_var}
#' \item \code{k_used}
#' \item \code{k_f}
#' \item \code{f}
#' \item \code{saturated}
#' \item \code{var_explained}
#' \item \code{PC}
#' \item \code{Vt}
#' \item \code{reconstruction}
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
#' calculate_f( zeroOutSmall(clustering_vector, 30), f_E=0.7, X_variances= matrixStats::colVars(brca) )
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
#' @param X l 
#' @param X_variances  l
#' @param clustering_vector l
#' @param clustering_list l
#' @param min_cl_size l
#' @param f_E l
#' @param k l
#' @return A \code{list} with following elements:
#' \itemize{
#' \item \code{C_base} A \code{clustering_list} representing input clustering.
#' \item \code{f} Recalculated \code{f_E} giving the fraction of variance of proper clusters to explain (excluding "noise" clusters (\code{C_0}))
#' \item \code{C_base_blockPCA} A \code{k}-th order clustering based reconstruction of \code{X}, based on \code{C_base}. It is output of \code{blockwisa_PCA_reconstruction} funtion.
#' }
#' @examples
#' data(brca)
#' data(brca_clusters)
#' result<- initial_clusterNreconstruct(X= brca, X_variances=matrixStats::colVars(brca),
#'				        clustering_vector=brca_clusters)
#' lapply(result$C_base, length)
#' length(result$C_base_blockPCA)
#' @export

initial_clusterNreconstruct<- function(X, 
				       X_variances,
				       clustering_vector,
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
	       return( list(C_base=clustering_list,
			    f=f,
			    C_base_blockPCA=NA) ) }
		
        # result if everything went OK	
     list(C_base=C_base, f=f,
     C_base_blockPCA=blockwise_PCA_reconstruction(X,
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
#' @param hclustering_list A list of clusters and sub-clusters indexed by strings of integer tuples (like '1','2','1,1','1,2','2,1','2,2' etc.). Each cluster is a vector of integers dentoing nodes belonging to that cluster.
#' @param cluster_tuple A vector of indexes describing the indexing tuple or a string (c(1,2) or '1,2')
#' @param subdivision A grouping of indexes inside cluster indexed by \code{cluster_tuple} into subgroups in form of \code{clustering_list}.
#' @return An updated \code{hclustering_list} with additional elements denoting subclusters of cluster indexed by \code{cluster_tuple}, indexed by the tuple convention ( e.g. if \code{cluster_tuple='1,2'} then added subclusters from subdivision with clusters '1','2','3' are '1,2,1', '1,2,2', '1,2,3').
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

#' @param hclustering_list A list of clusters and sub-clusters indexed by strings of integer tuples (like '1','2','1,1','1,2','2,1','2,2' etc.). Each cluster is a vector of integers dentoing nodes belonging to that cluster.
#' @return An integer giving the number of levels of hierarchy in \code{hclustering_list}
#' @export

get_max_lvl<- function(hclustering_list){
if (length(hclustering_list)==0) stop(" hclustering_list must be nonempty") 
names(hclustering_list)-> all_tuples
 tuple_lengths<-unlist( lapply(all_tuples, function(tpl) length(p.(tpl)) ))
 max(tuple_lengths)
}

#' Extract \code{clustering_list} describing division into the clusters at specified lvl of hierarchy of \code{hclustering_list}
#' 
#' @param hclustering_list A list of clusters and sub-clusters indexed by strings of integer tuples (like '1','2','1,1','1,2','2,1','2,2' etc.). Each cluster is a vector of integers denoting nodes belonging to that cluster.
#' @param g An integer specifying lvl of hierarchy (number of subdivision step) to extract from \code{hclustering_list}
#' @return A \code{clustering_list} describing division into the clusters at lvl \code{g} of hierarchy of \code{hclustering_list} 
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
#' @export

normalized_cl_entropy<- function( clustering_list ) {
 unlist(lapply(clustering_list, length))-> n_c
 p_c<- n_c/sum(n_c)
 length(n_c)-> n_groups
 -sum( p_c * log( p_c ) )/(log(n_groups))
}

#' Subdivide clusters in the innermost layer of hierarchical clustering and reconstruct the residuals in the subdivided clusters
#' 
#' @param X  l
#' @param X_similarity_matrix  l
#' @param X_variances  l
#' @param hclustering_list  l
#' @param hcluster_PCA_blocks l
#' @param k  l
#' @param f l
#' @param clfun2 l
#' @param clfun2OtherArgs_constant l
#' @param clfun2OtherArgs_ranges l
#' @param S_power l
#' @param corfun l
#' @return A \code{list}
#' @export


subclusterNreconstruct<- function(X, X_similarity_matrix, 
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
		S_Eg.C_k<-  abs (corfun( Eg [, C_K ])	)^S_power
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
#' @param X l
#' @param hcluster_blockPCA l
#' @param hclustering_list l
#' @param add_noise l
#' @param uncenter l
#' @return A matrix \code{Xg} which is the \code{g}-th order hierarchical clustering reconstruction of correlated variables in \code{X}, where \code{g} is the number of levels in \code{hclustering_list} and \code{hcluster_blockPCA}.

get_full_reconstruction<- function(X,hcluster_blockPCA, hclustering_list, add_noise=TRUE, uncenter=TRUE) {

g<-get_max_g(hclustering_list)
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


