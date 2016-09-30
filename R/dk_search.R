#' d, K Search Algorithm
#'
#' \code{dk_search} selects the optimal eigen-Laplacian dimensionality (\code{d}) and the number of clusters (\code{K})
#'   jointly based on model assessment criteria like Bayesian Information Criterion (BIC) and adjusted adherence.
#'
#' Details of the procedure is available in our manuscript.
#'
#' @references Qian Li, Damla Senturk, Catherine A. Sugar, Shanali Jeste, Charlotte DiStefano, Joel Frohlich, Donatello Telesca
#'   "\emph{Inferring Brain Signals Synchronicity from a Sample of EEG Readings}".
#' @param X_array list of data arrays, each of which is organized as No.objects * No.observations * No.segments.
#' @param max_d integer, maximal value of \code{d} and \code{K} to be considered
#' @param par.spectrum vector, spectral estimation parameters in the order of \code{c(a, wn, nn)} see \code{\link{spec.parzen}}
#' @param par.win vector, epoch smoothing parameters see \code{\link{ep_dismat}}
#' @param n.iter integer, number of iterations in \code{\link{MIC}} fitting
#' @return A list of objects with the following components:
#'   \item{\code{d}}{selected d}
#'   \item{\code{K}}{selected K}
#'   \item{\code{path}}{searching trajectory}
#' @examples
#' \dontrun{
#' # An example here
#'
#' ## Time series simulation:
#' ts_sim <- MIC_sim(alpha = 0.9, nsub = 3, segs = 10, fs = 100)$Data
#'
#'
#'
#' ## d,K searching: \strong{(approx. 10 mins running time)}
#' dk_search(ts_sim, n.iter = 2000, par.spectrum = c(50, 50), par.win = c(3,1))
#'
#' }
#'
#' @export
dk_search <- function(X_array,
                      max_d = 10,
                      n.iter,
                      par.spectrum = c(50),
                      par.win = c(1, 0)){
  result <- c()
  # Initiate
  D <- 2
  K <- 2
  while (D <= max_d){
    # data preparation for the new D
    list_data <- list()
    for (i in 1:length(X_array)){
      x <- X_array[[i]]
      list_data[[i]] <- MIC_prep(X = x, d = D, par.spectrum = par.spectrum, par.win = par.win)
    }
    if (D > 2) K <- K - 1 # conservative step-back
    output <- MIC(X = list_data, K = K, NumRun = n.iter)
    current_BIC <- output$ICs[[1]]
    current_COH <- output$ICs[[5]]
    result <- rbind(result, c(D, K, current_BIC, current_COH)) # track the trajectory

    while(K < max_d){
      output <- MIC(list_data, K = K + 1, NumRun = n.iter)
      result <- rbind(result, c(D, K + 1, output$ICs[[1]], output$ICs[[5]]))
      if (output$ICs[[1]]> current_BIC){
        current_BIC <- output$ICs[[1]]
        current_COH <- output$ICs[[5]]
        K <- K + 1
      } else {
        cat('D = ',D,': best_K = ',K,' COH = ',current_COH,'\n', sep = '')
        cat('-----------------------------------------','\n')
        break
      }
    }
    if(K == max_d){
      stop('Set a larger max_d!')
    }
    if (D >= K) break; # converge
    D <- K # next D
  }
  if(D != K){
    # Then chose D,K by max coherence, locally
    D_cand <- result[, 1] >= (D - 1)
    cand_result <- result[D_cand, ]
    pick <- which.max(cand_result[, 4])
    D <- cand_result[pick,1]
    K <- cand_result[pick,2]
  }
  return(list(d = D, K = K, path = result))
}
