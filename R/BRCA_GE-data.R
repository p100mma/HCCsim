#' BRCA gene expression dataset
#'
#' A subset of 500 random samples and 2500 random genes from BRCA dataset retrieved from European Genome-Phenome Archive, accesion ID: EGAS00000000083.
#'
#' @docType data
#'
#' @usage data(brca)
#'
#' @format A large numeric matrix where columns are genes and rows are samples.
#'
#' @references Pereira B, Chin   SF, Rueda   OM  et al.   The somatic mutation profiles of 2,433 breast cancers refine their genomic and transcriptomic landscapes. Nat Commun  2016;7:11479
#'
#' @source \url{https://ega-archive.org/studies/EGAS00000000083}
#' @examples 
#' data(brca)
#' dim(brca)
"brca"
