#' Cluster Aligner
#'
#' \code{clust_align} aligns two cluster assignments by relabelling the membership for
#'   a maximized concordance
#'
#' This aligner can take both \emph{vector} and \emph{matrix} clustering assignments, and aligns
#'   the second assignment maximally against the first assignment. It aims at reshuffling of
#'   cluster labels to facilitate a fair comparision between two clustering result.
#'
#' @param Z1 Clustering assignment No.1 (reference).
#' @param Z2 Clustering assignment No.2 (to be adjusted).
#' @param type Option of '\emph{vec}' and '\emph{mat}' as input format.
#' @return \code{Z2} aligned against \code{Z1} achieving maximal concordance.
#' @examples
#' Z1 <- c(1:5)
#' Z2 <- c(5:1)
#'
#' clust_align(Z1, Z2, type = 'vec')
#'
clust_align <- function(Z1, Z2, type = 'vec'){
  if (type == 'vec'){
    for (k in 1:length(unique(Z1))){
      Max <- sum(Z1 == k & Z2 == k) / (.01 + sum(Z2 == k) + sum(Z1 == k));
      for (tempk in  1:length(unique(Z2))){
        if (sum(Z1 == k & Z2 == tempk) / (.01 + sum(Z2 == tempk) + sum(Z1 == k)) > Max){
          Max <- sum(Z1 == k & Z2 == tempk) / (.01 + sum(Z2 == tempk) + sum(Z1 == k))
          dummy <- Z2 == k
          Z2[Z2 == tempk] <- k
          Z2[dummy] <- tempk
        }}}
    }
  if (type == 'mat'){
    for (k in 1:dim(Z1) [2]){
      for (tempk in  1:dim(Z2) [2]){
        Max <- sum(Z1 == Z2)
        Z2dummy <- Z2
        Z2dummy[, k] = Z2[, tempk]
        Z2dummy[, tempk] = Z2[, k]
        if (sum(Z1 == Z2dummy) > Max) Z2 <- Z2dummy
      }}}
  return(Z2)
}

#' Hard Clustering from MCMC samples
#'
#' \code{HardCluster} produces a point clustering estimate from a sample of posterior estimates (MCMC outputs)
#'
#' This procedure treats each clustering assignments as affinity matrices, and seek for a point estimate
#'   that minimizes the Frobenius norm to mean affinity matrix, from the assignments being traversed by the
#'   sample.
#'
#' @param ClustList A matrix of clustering results whose row vectors are assignments, that have been properly aligned.
#'   Usually organized as No.Iterations * No.Objects.
#' @return A list of objects with the following components:
#'   \item{\code{Cbest}}{Hard Clustering point estimate}
#'   \item{\code{vec}}{A vector of Frobenius norm to the mean affinity}
#'
#' @seealso \code{\link{MIC}} for its usage in MCMC, \code{\link{clust_align}} for clusters alignment.
HardCluster <- function(ClustList){
  numit <- dim(ClustList) [1]
  Ckern <- 0
  for (w in 1:numit){
    vec <- as.vector(ClustList[w, ])
    Ckern <- Ckern + outer(vec, vec, FUN = "==")
  }
  Ckern <- Ckern/numit
  countMax <- dim(ClustList)[2]^2+1
  dvec <- c()
  for (w in 1:numit){
    vec = as.vector(ClustList[w, ])
    countC = norm(Ckern - outer(vec, vec, FUN = '=='), 'F')
    if (countC < countMax){
      Cbest <- vec
      countMax <- countC
    }
    dvec <- append(dvec, countC)
  }
  list(Cbest = Cbest, vec = dvec)
}