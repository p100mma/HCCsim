
#' Compute squared error between two vectors
#'
#' @param x One of vectors between which squared error is to be computed. 
#' @param y One of vectors between which squared error is to be computed. 
#' @return A vector of squared errors between \code{x,y}.
#' @examples
#' SE( c(1,2,3), c(3,2,1) )
#' @seealso [MSE()]
#' @export

SE<- function(x,y) { (x-y)^2 }

#' Compute mean squared error between two vectors
#'
#' @param x One of vectors between which squared error is to be computed. 
#' @param y One of vectors between which squared error is to be computed. 
#' @return A mean of squared errors between \code{x,y}.
#' @examples
#' MSE( c(1,2,3), c(3,2,1) )
#' @seealso [SE()]
#' @export

MSE<- function(x,y) { mean(SE(x,y)) }

#' Compute weighted version of clustering coefficient from WGCNA utilizing matrix product
#'
#' This function produces clustering coefficient from similarity matrix just like \code{clusterCoef} from \code{WGCNA} package [2]. It is however (experimentally) much more time efficient because it explicitly involves a (matrix-matrix) matrix product. Speed up is drastic if R uses efficient matrix multiplication backend like OpenBLAS.
#'
#' @param A symmetric matrix of similarities between nodes, \code{A[i,j]} is weight of the connection between nodes \code{i,j}. All entries must be between \code{0} and \code{1}.
#' @return A vector of clustering coefficients - weighted version of the clustering coefficient designed for weighted correlation networks, as proposed by \code{WGCNA} methodology.
#' @examples
#' data(brca)
#' S<- similarity_matrix(fastPearsonData(brca))
#' cc<- ccfast(S)
#' @seealso \href{https://rmflight.github.io/posts/2023-11-01-installing-openblas-for-selfcompiled-r}{this link} [1] for how to set up OpenBLAS on Linux, 
#'	and \href{https://github.com/david-cortes/R-openblas-in-windows}{that one} for instructions for Windows.
#' @references [1]  M Flight, Robert. 2023. “Installing OpenBLAS for Self-Compiled R.” November 1, 2023. <https://rmflight.github.io/posts/2023-11-01-installing-openblas-for-selfcompiled-r>
#' @references [2] Langfelder P, Horvath S (2008). “WGCNA: an R package for weighted correlation network analysis.” BMC Bioinformatics, 559. <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559>
#' @export

ccfast<- function(A)    #optimized WGCNA clustering coefficient
      {
      diag(A)=0
      N<-nrow(A)
      A.A<-A %*% A
      numer<- vector()
      for(i in 1:N)
          numer[[i]]<- A[i,] %*% A.A[i,]
      denom= colSums(A)^2 - colSums(A^2)
      ifelse(denom==0,0,numer/denom)
      }

#' Compute strengths of nodes from similarity matrix
#'
#' @param A A similarity matrix  \code{A[i,j]} is weight of the connection between nodes \code{i,j}.
#' @return A vector of strengths of all nodes, weighted versions of degrees.
#' @export

strengths<- function(A) { diag(A)=0; colSums(A) }

#' Compute distance between empricial distribution functions
#'
#' Given and input value vector \code{y} and prepared \code{ecdf} object (empirical distribution function of other distribution), function will compute distance between empirical c.d.f. of \code{y} and the input \code{ecdf}.
#'
#' @param y A input vector of values which \code{ecdf} will be compared with the reference.
#' @param ecdf_x A function acting on \code{y}, reference empirical c.d.f. (or theoretical of a known distribution), e.g. produced by base R \code{ecdf()} function.
#' @param summarizeFun Optionally, an aggregation function called on the vector of absolute differences between two c.d.f.s.
#' @return If \code{summarizeFun} is \code{NULL} (default), a vector of absoulute differences between two compared c.d.f's computed on input \code{y}. If \code{summarizeFun} is given, return value of \code{summarizeFun} called on the mentioned vector of differences  is returned instead.
#' @examples
#' x<- rnorm(1000)
#' y<- rnorm(1000,0.05)
#' diffs<- ecdf_distance(x, ecdf(y))
#' mean_dist<- ecdf_distance(x, ecdf(y), mean)
#' @seealso [KS_ecdf_distance()] for specific version of this function
#' @export

ecdf_distance<- function(y, ecdf_x, summarizeFun=NULL) {
ecdf_y<- ecdf(y)
 distt=abs( ecdf_x(y) - ecdf_y(y) )
 if (!is.null(summarizeFun) ) distt=summarizeFun( distt) 
return(distt)
}

#' Kolmogorov-Smirnov distance between sample and reference c.d.f.
#'
#' @param y A input vector of values which \code{ecdf} will be compared with the reference.
#' @param ecdf_x A function acting on \code{y}, reference empirical c.d.f. (or theoretical of a known distribution), e.g. produced by base R \code{ecdf()} function.
#' @return A KS statistic computed in the following manner: \code{max(abs( ecdf_y(y) - ecdf_x(y)))}.
#' @examples
#' x<- rnorm(1000)
#' y<- rnorm(1000,0.05)
#' KS_dist<- KS_ecdf_distance(x, ecdf(y))
#' @seealso [ecdf_distance()] for more general distance between \code{ecdf}s.
#' @export

KS_ecdf_distance<- function(y, ecdf_x){
ecdf_distance(y,ecdf_x, max)
}
