#' Epoch-wise pairwise Distance Matrices
#'
#' \code{ep_dismat} calculates strucally smoothed pairwise distances across all objects.
#'   It is an essential prerequisite step of MIC implementation.
#'
#' This procedure consists of distance characterization across all objects to be clustered on
#'   and further smoothing into epochs.
#'
#' The distance metric is based on the \code{\link{TVD}} between a pair of spectral
#'   densities, calculated by the Parzen lag window estimator \code{\link{spec.parzen}}.
#'
#' Epochs are constructed by combining adjacent units (segments), which are further allowed
#'   to overlap with each other at certain amount.
#'
#' @param X A 3-D array of input organizing the data in the stucture of:
#'   No.Objects * No.Observations * No.Segments. Here objects refer to the clustering units;
#'   observations refer to the time series at each segment; and segments are pieces of time series
#'   that are repeatedly observed/measured of same length.
#' @param sf Sampling frequency
#' @param par.spectrum Parameter for \code{\link{spec.parzen}} in the order of \code{a, wn, nn}.
#' @param window Epoch-smoothing window configuration, in the order of (window size, overlap size).
#'   For example, non-smoothing setting makes \code{window = c(1,0)}.
#' @return a list object with the following components:
#'   \item{\code{diss_array}}{Pairwise distance array, No.objects*No.objects*No.epochs}
#'   \item{\code{ave_spec}}{Epochwise average spectral estimates using \code{\link{spec.parzen}}}
#'   \item{\code{fw}}{Frequencies used for spectral estimate}
#' @examples
#' Later on
#'
#' @export
ep_dismat <- function(X,
                      sf,
                      par.spectrum = c(100, sf/2, 512),
                      window = c(1,0)){
  if (length(dim(X)) != 3) stop("X needs to be a 3-D array!") else {
    dim(X)[1] -> nc
    dim(X)[2] -> nt
    dim(X)[3] -> ns
  }
  dt<-1/sf
  length.w<-512
  np<-length(par.spectrum)

  if (np==1) {a <- par.spectrum[1]; wn <- sf/2; length.w <- 512}
  if (np==2) {a <- par.spectrum[1]; wn <- par.spectrum[2]; length.w <- 512}
  if (np==3) {a <- par.spectrum[1]; wn <- par.spectrum[2]; length.w <- par.spectrum[3]}

  ## seg-wise normalized spectrum
  Spec <- apply(X, c(1,3), function(channel)
    spec.parzen(channel, a = a, dt = dt, wn = wn, nn = length.w)[, 2] / var(channel))
  fw <- spec.parzen(X[1,,1], a = a, dt = dt, wn = wn, nn = length.w)[, 1]

  ## Sliding window averaging.
  if(window[1]==1){
    a.Spec <- Spec
    ne <- ns
    incre = 1} else {
      incre <- window[1] - window[2]
      ne <- floor((ns-window[1])/incre)+1
      a.Spec <- array(NA, dim = c(length.w, nc, ne))
      for(k in 1:ne){
        start <- (k - 1) * incre + 1
        end <- start + window[1] - 1
        a.Spec[, , k]<-apply(Spec[, , start:end], c(1, 2), mean)
        }
    }
  rm(Spec)

  ## Epoch distance mat
  MatDiss<-array(NA, c(nc, nc, ne))
  for (k in 1:ne){
    for (i in 1:nc) for (j in i:nc) {
      MatDiss[i, j, k] <- TVD(w = fw, a.Spec[, i, k], a.Spec[, j, k])
      MatDiss[j, i, k] <- MatDiss[i, j, k]
      }
    diag(MatDiss[, , k]) <- 0
    Mat <- MatDiss[,,k]
    Mat [lower.tri(Mat)] <- t(Mat) [lower.tri(Mat)] #Symmetry
    MatDiss[, , k] <- Mat
    rm(Mat)
  }
  return(list(diss_array = MatDiss, ave_spec = a.Spec, fw = fw))
}