#' Eigen-Laplacian Dimension Reduction
#'
#' \code{eigen_lap} carries out eigen-decomposition on the Laplacian matrices, which
#'   operates on distance matrices and yields low dimensional representations of original objects.
#'
#' This procedure is a direct implementation of the spectral clustering technique explicated by [Ng 2001].
#' It is also a preprocessing step for \code{\link{MIC}}, following distance characterization by \code{\link{ep_dismat}}.
#'
#' @references Ng, Andrew Y., Michael I. Jordan, and Yair Weiss. \emph{On spectral clustering: Analysis and an algorithm.}
#'   Advances in neural information processing systems 2 (2002): 849-856.
#'
#' @param X Output \code{diss_array} from \code{\link{ep_dismat}}, featuring pairwise distances on epoch-smoothed units.
#' @param D Dimensionality of the output, which correspond to the number of eigenvectors extracted from graph Laplacian.
#' @return A list of objects with the following components:
#'   \item{\code{eig_data}}{A 3-D array of the Eigen-Laplacian representation.}
#'   \item{\code{eig_value}}{Corresponding eigen-values of the extracted eigen-vectors}
#' @examples
#' Later on
#'
eigen_lap <- function(X, D){
  nc <- dim(X) [1]
  ne <- dim(X) [3]
  eig_data <- array(NA, c(nc, D, ne))
  eig_val <- matrix(NA, nrow = ne, ncol = nc)
  for (j in 1:ne){
    simMat <- 1 - X[, , j]
    diag(simMat) <- 0
    dia_mat <- diag(1 / sqrt(rowSums(simMat)))
    L <- dia_mat %*% simMat %*% dia_mat
    L_eig <- eigen(L, symmetric = TRUE)
    #vec <- L_eig$vectors
    eig_val [j, ] <- L_eig$values
    eigscore <- L_eig$vectors [, 1:D]
    Ln <- diag(eigscore %*% t(eigscore))
    eig_data[, , j] <- diag(1 / sqrt(Ln)) %*% eigscore
    rm(simMat, dia_mat, L, L_eig, eigscore, Ln)
  }
  return(list(eig_data = eig_data, eig_value = eig_val))
}