# Helper Functions in MIC package

#' Total Variation Distance
#'
#' \code{TVD} returns the Total Variation Distance between two densities.
#'
#' @param w A vector of sampled points.
#' @param f1 Density 1 evaluated at values \code{w}.
#' @param f2 Density 2 evaluated at values \code{w}.
#' @return A scalar, in range [0,1] if both \code{f1} and \code{f2} are normalized densities.
#' @examples
#' w <- seq(0, 1, by =.05)
#' f1 <- dbeta(w, 2, 2); f2 <- dbeta(w, 2, 5)
#'
#' TVD(w, f1, f2)
#' TVD(w, f2, f1)
#'
TVD <- function(w, f1, f2){
  #Compute TV distance using the trapezoidal rule
  if (length(w) != length(f1) | length(w) != length(f2)) stop("w, f1 and f2 must have the same length")
  n <- length(w)
  int <- (w[2:n]-w[1:(n-1)])%*%(pmin(f1,f2)[2:n]+pmin(f1,f2)[1:(n-1)])/2
  return(as.double(1-int))
}

#' Estimate Spectral Density using Parzen Lag Window
#'
#' \code{spec.parzen} calculate the spectral estimate based on a Fourier transform of a
#'   truncated and Parzen window smoothed auto-covariance function (ACF).
#'
#' The raw periodogram \code{\link{spec.pgram}} is not a consistent estimator of the spectral density,
#'   therefore a class of lag window estimators are considered as surrogates in practice achieving
#'   smootheness and consistency.
#'
#' Parzen window estimator works the best when the true spectrum is continuous, and specifically
#'   have peak concentrations at certain frequencies. Such estimators operates on times series
#'   that are presumably zero-mean and stationary, so that demean and detrending are highly
#'   recommended on \code{x} before implementation.
#'
#' @param x A univariate time series.
#' @param a Max lag to truncate the sample ACF estimates, default \code{a=100}. It cannot exceed
#'   the number of observations in \code{x}.
#' @param dt 1/freq, unit sampling interval of \code{x}.
#' @param w0 The minimal frequency of interest.
#' @param wn The maximal frequency of interest.
#' @param nn Resolution of the spectral estimates (number of points on frequency domain).
#' @return A \code{nn}-by-2 matrix, with the first column being the estimated frequencies,
#'   and the second column being the esimated spectrum.
#' @examples
#' x <- rnorm(100)
#' spec <- spec.parzen(x, a = 50, nn = 50)
#'
#' plot(spec, type = 'l')


spec.parzen <- function( x,
                         a = 100,
                         dt = 1/length(x),
                         w0 = 10^(-5),
                         wn = 1/(2*dt),
                         nn = 512){
  dt <- dt*(2*pi)
  # Compute the smoothed periodogram using a Parzen window
  if (a >= (length(x) - 1)) stop("The bandwidth value is too big")
  ## parzen window smoother
  kp <- numeric(2 * a - 1)
  tt <- seq(-1, 1, length = 2*a - 1)
  kp[abs(tt) < .5] <- 1 - 6*abs(tt[abs(tt) < .5])^2 + 6*(abs(tt[abs(tt) < .5])^3)
  kp[abs(tt) >= .5]<-2*(1 - abs(tt[abs(tt) >= .5]))^3
  ## acf
  Cov <- as.vector(acf(x, lag.max = a-1, type = "covariance", plot=F)$acf)
  CovS <- kp[a: (2*a-1)] * Cov
  time <- seq(0, length(x) * dt, by = dt)[1: a]
  w <- seq(w0, wn, length.out = nn)
  SS <- sapply(w, function(x) (2 * dt / pi) * sum(cos(x * time)[-1] * CovS[-1]) + (dt / pi) * CovS[1])
  Result <- matrix(c(w,SS), ncol = 2)
  colnames(Result) <- c('freq', 'spec')
  return(Result)
}
