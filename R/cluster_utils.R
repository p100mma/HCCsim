
unlapply<-function(...) unlist(lapply(...))

#' A minimal, efficient table function for one dimensional vectors
#' 
#' Emulates \code{table} in an efficient manner for a narrow case of one dimensional atomic vectors. Note that no correctness checks are done.
#' 
#' @param values A 1 dimensional atomic vector
#' @return A \code{list} wirth the following elements: 
#' \itemize{
#' \item \code{value} a vector of unique values in \code{values}.
#' \item \code{count} a vector of counts of each element of \code{value} component, i.e. \code{fastTable(x)$count[i]} is frequency of \code{fastTable(x)$value[i]} in \code{x}.
#' 	   }
#' @examples
#' fastTable( sample(c(1,2,3,4),100,replace=TRUE,prob=c(0.2,0.1,0.5,0.2)))
#' @export
fastTable<- function(values)
{
    uqs<-unique(values)
    histo<-unlist(lapply(uqs, function(x) sum(values==x)))
    list(value=uqs,
         count=histo)
}


#' Relabel small groups to \code{0}
#'
#' Based on input label vector, zero out labels corresponding to groups smaller than given minimal group size.
#'
#' @param labelV An integer membership vector for which \code{cl_mem[i]} gives a label of the group to which node \code{i} belongs. A value of \code{0} by assumption encodes nodes that will be not taken into the account by the relabeling operation.
#' @param mS Integer scalar giving minimum group size that is not relabeled to \code{0}.
#' @return An altered integer membership vector in which each label corresponding to the group of size \code{ < mS } is changed to \code{0}.
#' @examples
#' mem<- c(0,0,0,1,0,1,2,3,0,1,2,2,2,1,4,4,0,0)
#' zeroOutSmall(mem,3)
#' @export
zeroOutSmall<- function(labelV, mS) {
rment<- labelV
non0labs<- unique(labelV[labelV!=0])
non0cnts<- unlapply(non0labs, function(lab) sum(labelV==lab))
for (l in seq_along(non0labs)) 
    if (non0cnts[[l]] < mS )
       rment[ labelV==non0labs[[l]] ]=0
return(rment)
}


#' Relabel group labels to consecutive integers
#'
#' Based on input integer label vector, normalize labels such that if there are \code{M} unique labels in total, each of them gets mapped to one of the integers in \code{1:M}.
#'
#' @param labelV An integer membership vector for which \code{cl_mem[i]} gives a label of the group to which node \code{i} belongs. A value of \code{0} by assumption encodes nodes that will be not taken into the account by the relabeling operation.
#' @return An altered integer membership vector in which each label corresponding to the group of size \code{ < mS } is changed to \code{0}.
#' @examples
#' mem<- c(0,0,0,1,0,1,67,3,0,1,67,67,5,5,4,4,0,0)
#' tidyUpLabels(mem)
#' @export
tidyUpLabels<- function(labelV) {
uqn0<-unique(labelV[labelV!=0])
rmentL<- seq_along(uqn0)
rment<- rep(0, length(labelV))
attributes(rment)<- attributes(labelV)
for (l in seq_along(uqn0))
   rment[ labelV == uqn0[[l]] ]= rmentL[[l]]
return(rment)
}

new_HCCSim_clustering_vector<- function(x, params=list() )  {
stopifnot(is.integer(x))
structure(x,
	  class="HCCSim_clustering_vector",
	  params= params
	 )

}

check_HCCSim_clustering_vector<- function(x) {
vec<- unclass(x)
if (any(is.na(vec))) stop("no NA labels allowed")
if (any(vec<0)) stop("labels inside clustering_vector can only be 0 or positive integers")
x
}

#' Initialize clustering_vector object
#'
#' \code{clustering_vector} is an integer vector for which \code{clustering_vector[i]} gives a label of the group to which node \code{i} belongs. 
#'  It posesses additional attribute \code{'params'} which stores optional additional parameters of the clustering algorithm used.
#'  Labels of proper clusters start from 1. Zero label denotes the "noise" subset of nodes not belonging to any proper cluster.
#' 
#' @param x A vector of nonnegative integer labels. Gets coerced to \code{integer} type.
#' @param params An arbitrary \code{list} that will be stored in an attribute of resulting \code{clustering_vector}.
#' @return A \code{HCCSim_clustering_vector} object.
#' @examples
#' data(brca_clusters)
#' brca_clusters<- HCCSim_clustering_vector(brca_clusters, params=list(algorithm='mcl', inflation=2))
#' @export
HCCSim_clustering_vector<- function(x= integer(), params=list()) {
x = as.integer(x)
check_HCCSim_clustering_vector(new_HCCSim_clustering_vector(x, params=params))
}


new_HCCSim_clustering_list<- function(x=list(),
				      domain_size=length(unlist(x)),
				      params=list() )  {
structure(x,
	  class="HCCSim_clustering_list",
	  domain_size=domain_size,
	  params= params
	 )
}

check_HCCSim_clustering_list<- function(x) {
unlisted<- unlist(x)
if ((length(unlisted) > attr(x,'domain_size'))) stop("total number of elements inside each cluster must be less or equal to domain_size")
if ((length(unlisted) > length(unique(unlisted)))) stop("Clusters must be non-overlapping, no nodes should repeat between clusters")
if (any(is.na(unlisted))) stop("no NA indexes allowed")
if (any(!(unlisted %in%  1:attr(x,'domain_size')))) stop("Nodes must have indexes in 1:domain_size")
x
}

#' Initialize clustering_list object
#'
#' \code{clustering_list} is a list of named elements. Each element represents one proper cluster of nodes and contains integers denoting nodes belonging to that cluster.
#' Name of the element is the cluster label (integer converted to string, e.g. for two proper clusters, there are two elements '1' and '2').
#' Note that info on clusterless, "noise" nodes is deduced by taking into the account which integers are not contained inside proper clusters.
#' This requires using attribute \code{'domain_size'} to store the total number of nodes in the network which was clustered.
#'  It posesses additional attribute \code{'params'} which stores optional additional parameters of the clustering algorithm used.
#' A more friendly function to initialize such object is \code{clustering_vector2list} which accepts as an input a more typical in R representation of clusters, membership vector and converts it to \code{clustering_list}.
#' 
#' @param x A list of clusters. Each element is a vector of indexes of nodes that belong to each cluster. Names denote cluster labels. Can exclude some nodes of originally clustered network, which makes them clusterless.
#' @param domain_size Total number of nodes in the network of which node indexes have been placed in each cluster.
#' @param params An arbitrary \code{list} that will be stored in an attribute of resulting \code{clustering_list}.
#' @return A \code{HCCSim_clustering_list} object.
#' @examples
#' data(brca_clusters)
#' labels<- sort( unique(brca_clusters[brca_clusters!=0]))
#' input_list<-list()
#' for (i in labels)
#'    input_list[[i]]<- which(brca_clusters==i)
#' HCCSim_clustering_list(x= input_list, domain_size=length(brca_clusters),  params=list(algorithm='mcl', inflation=2)) 
#' @export

HCCSim_clustering_list<- function(x=list(), domain_size=length(unlist(x)), params=list()) {
clist<- new_HCCSim_clustering_list(x=x,
				   domain_size=domain_size,
				   params=params)
check_HCCSim_clustering_list(clist)
}


#' Transform clustering_vector to clustering_list
#'
#' \code{clustering_list} is a list of named elements. Each element represents one proper cluster of nodes and contains integers denoting nodes belonging to that cluster.
#' Name of the element is the cluster label (integer converted to string, e.g. for two proper clusters, there are two elements '1' and '2').
#' Note that info on clusterless, "noise" nodes is deduced by taking into the account which integers are not contained inside proper clusters.
#' This requires using attribute \code{'domain_size'} to store the total number of nodes in the network which was clustered.
#'  It posesses additional attribute \code{'params'} which stores optional additional parameters of the clustering algorithm used.
#' \code{clustering_vector} is an integer vector for which \code{clustering_vector[i]} gives a label of the group to which node \code{i} belongs. 
#'  It posesses additional attribute \code{'params'} which stores optional additional parameters of the clustering algorithm used.
#'  Labels of proper clusters start from 1. Zero label denotes the "noise" subset of nodes not belonging to any proper cluster.
#'
#' @param clustering_vector An integer vector for which \code{clustering_vector[i]} gives a label of the group to which node \code{i} belongs. Its attribute \code{'params'} is copied onto the resulting \code{clustering_list}. Does not have to be an object of class \code{HCCSim_clustering_vector} but function assumes it sticks to the labeling convention outlined for that kind of representation.
#' @return A list with named elements with one element per proper cluster, where each element is integer vector of indexes of nodes. An object of class \code{HCCSim_clustering_list}.
#' @examples
#' data(brca_clusters)
#' clustering_vector2list( HCCSim_clustering_vector(brca_clusters, params= =list(algorithm='mcl', inflation=2)) )
#' @export


clustering_vector2list<- function(clustering_vector){

labels<- sort( unique(clustering_vector[clustering_vector!=0]))
init_list<-list()
for (i in labels)
	init_list[[as.character(i)]]<- which(clustering_vector==i)
HCCSim_clustering_list(init_list, domain_size= length(clustering_vector), params= attr(clustering_vector, 'params') )
}

#' Transform clustering_list to clustering_vector
#'
#' \code{clustering_list} is a list of named elements. Each element represents one proper cluster of nodes and contains integers denoting nodes belonging to that cluster.
#' Name of the element is the cluster label (integer converted to string, e.g. for two proper clusters, there are two elements '1' and '2').
#' Note that info on clusterless, "noise" nodes is deduced by taking into the account which integers are not contained inside proper clusters.
#' This requires using attribute \code{'domain_size'} to store the total number of nodes in the network which was clustered.
#'  It posesses additional attribute \code{'params'} which stores optional additional parameters of the clustering algorithm used.
#' \code{clustering_vector} is an integer vector for which \code{clustering_vector[i]} gives a label of the group to which node \code{i} belongs. 
#'  It posesses additional attribute \code{'params'} which stores optional additional parameters of the clustering algorithm used.
#'  Labels of proper clusters start from 1. Zero label denotes the "noise" subset of nodes not belonging to any proper cluster.
#'
#' @param clustering_list A list of named elements representing division of nodes into clusters, assumed to be built with the convention of \code{HCCSim_clustering_list}. Function does not require that argument to be object of the mentioned class but assumes it posseses the same attributes as such object and is structured similarly (attribute \code{domain_size} being crucial!).
#' @return An integer vector for which \code{clustering_vector[i]} gives a label of the group to which node \code{i} belongs. An object of class \code{HCCSim_clustering_vector} with attribute \code{params} copied from the input \code{clustering_list}.
#' @examples
#' data(brca_clusters)
#' cl_vec<-HCCSim_clustering_vector(brca_clusters, params= =list(algorithm='mcl', inflation=2)) 
#' table(cl_vec)
#' clustering_vector2list( cl_vec )-> cl_list
#' lapply(cl_list, length)
#' clustering_list2vector(cl_list) -> cl_vec2
#' all( cl_vec == cl_vec2)
#' @export
clustering_list2vector<- function(clustering_list){
init_vector<- rep(0, attr(clustering_list, 'domain_size'))
for (i in seq_along(clustering_list))
	init_vector[ clustering_list[[i]] ] = i
HCCSim_clustering_vector(init_vector, params= attr(clustering_list, 'params'))
}
