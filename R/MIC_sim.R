#' Time Series Simulation for MIC
#'
#' \code{MIC_sim} simulates multivariate picewise stationary time series for a sample of subjects.
#'
#' @references Qian Li, Damla Senturk, Catherine A. Sugar, Shanali Jeste, Charlotte DiStefano, Joel Frohlich, Donatello Telesca
#'   "\emph{Inferring Brain Signals Synchronicity from a Sample of EEG Readings}".
#'
#' @param alpha numeric in [0, 1], group-subject adherence
#' @param nsub integer, number of subjects
#' @param segs integer, number of segments observed
#' @param fs integer, sampling frequency
#' @param Ct vector of integers, true group labels
#' @param bad_sub integer, number of outlying subjects
#' @param scheme stationarity: 0 = stationary; 1 = piecewise stationary
#'
#' @return A list of objects with the following components:
#'   \item{\code{C}}{True group label}
#'   \item{\code{Ci}}{True individual label}
#'   \item{\code{Data}}{Simulated Data, list of 3-D arrays}
#'   \item{\code{Points}}{Transition time points}
#'
#' @examples
#' # Stationary simulation:
#'
#' x_stat <- MIC_sim(alpha = 0.9, nsub = 2, segs = 10, fs = 100)
#'
#' # # sample adherence (alpha)
#'
#'   sum(x_stat$C == x_stat$Ci[, 1]) / length(x_stat$C)
#'
#' # # stationarity
#'
#'   x_stat$Points
#'
#' # Non-stationary (piecewise stationary) simulation:
#'
#' x_nstat <- MIC_sim(alpha = 0.9, nsub = 2, segs = 10, fs = 100)
#' x_nstat$Points
#'
#' @export
MIC_sim <- function(alpha,
                    nsub,
                    segs,
                    fs, Ct = rep(c(1, 2, 3, 4), rep(10, 4)), bad_sub = 0, scheme = 0){
  # Helper function
  pars<- function(eta, M = 1.1, fs){
    phi1 <- - 1 / M ^ 2
    phi2 <- 2 * cos(2 * pi * eta / fs) / M
    return(c(phi2, phi1))
  }
  index <- c(1:length(unique(Ct)))
  nelec <- length(Ct)
  ## Subject level labels
  sub_label <- matrix(NA, nrow = nelec, ncol = nsub)
  for (i in 1:nsub){
    sub_label[, i] <- Ct
    if (i <= bad_sub){
      change <- runif(length(Ct)) > 0.5}
    else{
      change <- runif(length(Ct)) > alpha}
    for (j in 1:length(Ct)){
      if (change[j]){
        sub_label[j, i] <- index[ -Ct[j]][sample(length(index) - 1, 1)]
      }
    }
  }
  data <- list()

  ## Non-stationary scheme:
  ##  0: stationary
  ##  1: exponential + normal
  breaks <- list()
  for (i in 1:nsub){
    if(scheme == 0){
      breaks[[i]] <- c(0, segs * fs)
      } else if (scheme == 1){
      breaks[[i]] <- c(0)
      step = 0
      while(tail(breaks[[i]], 1) < (segs * fs)){
        if(step == 0){
          last = tail(breaks[[i]], 1)
          breaks[[i]] <- c(breaks[[i]],last + floor(rexp(1, 0.05) * fs))
          step <- 1
          next
        }
        else{
          last = tail(breaks[[i]], 1)
          breaks[[i]] <- c(breaks[[i]], last + floor(rnorm(1,mean = 5,sd = 1) * fs))
          step <- 0
        }
      }
      breaks[[i]] <- c(breaks[[i]][-length(breaks[[i]])], segs * fs)
    } else {stop('No such scheme!')}
  }
  weights1 <- cbind(c(1,2,0,0,0),c(0,1,2,0,0),c(0,0,1,1,0),c(0,0,0,1,1))
  weights2 <- cbind(c(1,0,0,2,0),c(0,0,2,1,0),c(0,0,1,2,0),c(0,0,0,1,2))
  for(sub in 1:nsub){
    mat_data <- matrix(nrow = nelec, ncol = segs * fs)
    before <- rep(rep(c(T, F), length.out = length(breaks[[i]]) - 1), diff(breaks[[i]]))
    after <- !before
    for(e in 1:nelec){
      ts1 <- arima.sim(n = segs*fs, model = list(ar = pars(2, fs = fs)))
      ts2 <- arima.sim(n = segs*fs, model = list(ar = pars(6, fs = fs)))
      ts3 <- arima.sim(n = segs*fs, model = list(ar = pars(10, fs = fs)))
      ts4 <- arima.sim(n = segs*fs, model = list(ar = pars(21, fs = fs)))
      ts5 <- arima.sim(n = segs*fs, model = list(ar = pars(40, fs = fs)))
      ts <- rbind(ts1,ts2,ts3,ts4,ts5)
      mat_data[e,before] <- t(weights1[,sub_label[e,sub]]) %*% ts[,before]
      mat_data[e,after] <- t(weights2[,sub_label[e,sub]]) %*% ts[,after]
      mat_data[e,] <- mat_data[e,]
    }
    dim(mat_data) <- c(nelec, fs, segs)
    data[[sub]] <- mat_data
  }
  return(list(C = Ct, Ci = sub_label, Data = data, Points = breaks))
}