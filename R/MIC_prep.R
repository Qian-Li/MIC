#' Data preparation for MIC
#'
#' \code{MIC_prep} prepares stuctured data for \code{MIC}. It sequentially carries out signal detrending, epoch smoothing,
#'   distance characterization and dimensional reduction sequentiall.
#'
#' @param X 3-d array, organized as No.objects * No.observations * No.segments, see simulated example by \code{\link{MIC_sim}}
#' @param d integer, dimensionality of the eigen-Laplacian representation
#' @param exclude vector of integers, indicating objects to be excluded
#' @param par.spectrum vector, spectral estimation parameters in the order of \code{c(a, wn, nn)} see \code{\link{spec.parzen}}
#' @param par.win Vector, epoch smoothing parameters see \code{\link{ep_dismat}}
#' @return List of data matrices, each with No.objects rows and \code{d} columns.
#' @examples
#' \dontrun{
#' # Simulated data:
#' ts_sim <- MIC_sim(alpha = 0.9, nsub = 3, segs = 10, fs = 100)$Data
#'
#' # Data preparation, subject 1 epoch 1
#' sub1 <- MIC_prep(ts_sim[[1]], d = 3, par.spectrum = c(50, 50), par.win = c(3, 1))
#'
#'
#' # No. of epochs
#' length(sub1)
#'
#'
#' # Epoch data: No.objects * d
#' dim(sub1[[1]])
#' }
#'@seealso \code{\link{MIC}} for its usage and \code{\link{MIC_sim}} for time series simulation.
#'
#' @export
MIC_prep <- function(X, d,
                     exclude = c(),
                     par.spectrum = c(50),
                     par.win = c(1, 0)){
  if (length(exclude) != 0) X <- X[- exclude, , ] #exclusion
  segs <- dim(X) [3]
  fs <- dim(X) [2]
  nc <- dim(X) [1]
  # Detrending
  for (i in 1:nc){
    for (j in 1:segs){
      lmodel <- lm(X[i, , j] ~ c(1:fs))
      X[i, , j] <- X[i, , j] - predict(lmodel)
    }
  }
  diss_out <- ep_dismat(X, sf = fs, par.spectrum = par.spectrum, window = par.win)$diss_array
  spec_out <- eigen_lap(diss_out, d)$eig_data
  list_out <- lapply(seq(dim(spec_out)[3]), function(xx) spec_out[, , xx])
  rm(diss_out, spec_out)
  return(list_out)
}