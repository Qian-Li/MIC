#' Cluster Aligner
#'
#' \code{clust_align} aligns two cluster assignments by a maximally concordant relabeling
#'
#' This aligner can take both \emph{vector} and \emph{matrix} clustering assignments, and aligns
#'   the second assignment against the first assignment. It reshuffles the second assignment labels
#'   such that two assignemnts are maximally matched. This procedure facilitates a direct and
#'   fair comparision between two clustering result.
#'
#' @param Z1 vector or matrix, clustering assignment 1 (reference).
#' @param Z2 vector or matrix, clustering assignment 2 (to be relabeled).
#' @param type what type of input format: "\code{vec}" or "\code{mat}"
#' @return \code{Z2}, aligned against \code{Z1} with maximal concordance.
#' @examples
#' Z1 <- c(1:5)
#' Z2 <- c(5:1)
#'
#' clust_align(Z1, Z2, type = 'vec')
#'
#' @export

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
#' \code{HardCluster} produces a point clustering estimate from a sample clustering labels (MCMC outputs)
#'
#' This procedure treats each clustering assignments as affinity matrices, and seek for a point estimate
#'   that minimizes the Frobenius norm to mean affinity matrix, from the assignments being traversed by the
#'   sample.
#'
#' @param ClustList matrix, with row vectors being cluster labels that are properly aligned.
#'   Usually organized as No.iterations * No.objects.
#' @return A list of objects with the following components:
#'   \item{\code{Cbest}}{vector, hard clustering estimate}
#'   \item{\code{vec}}{vector, Frobenius norm to the mean affinity (analog of deviation from the mean)}
#'
#' @seealso \code{\link{MIC}} for its usage in MCMC, \code{\link{clust_align}} for clusters alignment.
#'
#'
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

#' MIC Helper: logSum
#'
#' \code{logSum} is a samll helper function calculating the sum of logarithms.
#'
#' @param l vector
#' @return numeric
#' @seealso \code{\link{MIC}} for its usage
logSum <- function(l) {max(l) + log(sum(exp(l - max(l))))}

#' MIC Helper: margPi
#'
#' \code{margPi} is a samll helper function calculating the marginal prior of MIC.
#'
#' @param p vector of proportions, sum up to 1
#' @param a numeric, adherence
#' @return numeric
#' @seealso \code{\link{MIC}} for its usage

margPi <- function(p,a){return(p * a + (1 - p) * (1 - a) / (length(p) - 1))}

#' MIC_sim Helper: pars
#'
#' \code{pars} determines the AR(2) coefficient with desired oscillation properties: peak location, peak width and
#'  sampling frequency.
#'
#' @param eta numeric, peak location
#' @param M numeric greater than 1, narrower as M approximates 1
#' @param fs sampling frequency
#' @return vector of AR(2) coefficients
#' @seealso \code{\link{MIC_sim}} for its usage
#'
pars<- function(eta, M = 1.1, fs){
  phi1 <- - 1 / M ^ 2
  phi2 <- 2 * cos(2 * pi * eta / fs) / M
  return(c(phi2, phi1))
}