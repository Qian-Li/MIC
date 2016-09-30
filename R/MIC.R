#' Multilevel Integrative Clustering
#'
#' \code{MIC} implements integrative clustering on highly structured data. The admissible structure contains repeated measurements
#'   from multiple subjects, where both a group clustering and individual clusterings are of inferential interests.
#'
#' @param X list of lists, and each sublist contains data matrices of the same size, for a proper format see
#'   the output of \code{\link{MIC_prep}}
#' @param K integer, number of clusters
#' @param NumRun number of iterations for MCMC sampling, with 1/5 used for burn-in
#' @param a numeric, hyperparameter for adherence ~ TBeta(a, b)
#' @param b numeric, hyperparameter for adherence ~ TBdta(a, b)
#' @param IndivAlpha boolean, whether to assume individual adherence (population - subjects)
#' @param IndivBeta boolean, whether to assume individual adherence (subjects - epochs)
#' @param Concentration numeric, hyperparameter for Dirichlet prior
#' @param adj numeric, adjustment term for label sampling with default of \code{5e-05}
#' @param block integer, unit size of MCMC samples to be written externally (\strong{memory saving!})
#' @param time.start the start of code timing
#' @return A list of objects with the following components:
#'   \item{\code{Alpha, AlphaBounds}}{point estimates and CIs of subject adherence to a group estimate}
#'   \item{\code{Betas, BetaBounds}}{point estimates and CIs of epoch adherence to its subject estimate}
#'   \item{\code{Cbest}}{point estimate of group level clustering result}
#'   \item{\code{Cibest}}{point estiamte of subject level clustering result}
#'   \item{\code{Lbest}}{point estimates of epoch level clustering result}
#'   \item{\code{ICs}}{model assessment metrics: BIC, BIC2, DIC4, DIC7, adjusted coherence}
#'   \item{\code{...}}{variability measures, posterior clustering probability matrices, etc. Please see our manuscript for details}
#' @references Qian Li, Damla Senturk, Catherine A. Sugar, Shanali Jeste, Charlotte DiStefano, Joel Frohlich, Donatello Telesca
#'   "\emph{Inferring Brain Signals Synchronicity from a Sample of EEG Readings}".
#' @seealso \code{\link{eigen_lap}}, \code{\link{ep_dismat}} and \code{\link{MIC_prep}} for data preparation, \code{\link{spec.parzen}} for the estimation
#'   of spectrum and \code{\link{TVD}} for pairwise distance metric, \code{\link{HardCluster}} for point clustering estimate from
#'   MCMC samples.
#' @examples
#'
#' \dontrun{
#'   # Time Series simulation:
#'   ts_sim <- MIC_sim(alpha = 0.9, nsub = 3, segs = 10, fs = 100)
#'
#'
#'   # Data preparation:
#'   list_data <- list()
#'   for (i in 1:3) list_data[[i]] <- MIC_prep(X = ts_sim$Data[[i]], d = 4,
#'     par.spectrum = c(50, 50, 100), par.win = c(3, 1))
#'
#'
#'
#'   # MIC: (Running time: 3mins)
#'   output <- MIC(X = list_data, K = 4, NumRun = 5000)
#'
#'
#'
#'   # Clustering accuracy: group level
#'   sum(clust_align(ts_sim$C, output$Cbest, type = 'vec') == ts_sim$C) / 40
#'   ## [1] 1
#'
#'
#'   # Clustering accuracy: individual level
#'   for (i in 1:3){
#'   aligned_label <- clust_align(ts_sim$Ci[, i], output$Cibest[, i], type = 'vec')
#'   cat(sum(aligned_label == ts_sim$Ci[, i]) / 40, '\n')
#'   }
#'   ## 1 1 1
#' }
#'
#' @export
MIC <- function(X,K,NumRun,# Required pars
                a = 1, b = 1,
                IndivAlpha = TRUE,
                IndivBeta = TRUE,
                Concentration = 1,
                adj = 5e-05,
                block = 1000,
                time.start = Sys.time()){
  # ------ MIC starts
  cat('--- New MIC session: ---','\n')
  # ------ Initialization
  #	# Hyper pars
  Gamma <- rep(1, K) #Dirichlet concentration
  a0 <- list() ## Optional input
  b0 <- list() ## Optional input
  mu0 <- list()## Optional input

  #	# Fixed pars
  N <- length(X) #Number of subjects
  M <- sapply(X, length) #Number of epochs by subject
  Total_e <- sum(M)
  P <- dim(X [[1]][[1]]) [1] #subjects total to be clustered on
  D <- rapply(X, ncol, how = "list") #Vector of dimension for each data source
  S <- list() #sample variance for each data source

  # # Output pars
  Lbest <- list() #List of best clustering for each data source
  Cibest <- matrix(NA, nrow = P, ncol = N) #List of best hard cluster for each batch
  ProbC <- matrix(0, nrow = P, ncol = K)
  ProbCi <- list()
  ProbL <- list()

  #	# Intermediate pars
  L <- list() #Current clustering for each data source (updates each MCMC draw)
  Ci <- list() #Current clustering for each batch source (updates each MCMC draw)
  Llist <- list() #Temp output of epoch clustering result
  Cilist <- list() #Temp output of subject clustering result
  Clist <- matrix(NA, nrow = P, ncol = block) #Temp output of population clustering result
  logC <- list() #log clustering probabilities for each data batch
  Lp <- list() #clustering probabilities for each data source
  A <- list() #Posterior gamma shape parameter
  B <- list() #Posterior gamma rate parameter
  #  Tau = list() #1/Tau = posterior variance, no need to record
  Sigma <- list() #Posterior variance
  mu <- list() #List of means for each data source, by cluster(update each draw)
  mpi <- matrix(nrow = N, ncol = K) # Mixing probability lvl1
  mpij <- list() # Mixing probaility lvl2
  Cprob <- matrix(nrow= P, ncol= K) #overall clustering probabilities

  Betalist <- vector("list", N) #List of beta(batch and source adherence) parameters, for cbind later.
  nu_e <- list()
  Csize <- list()
  beta <- list()

  # IC outputs
  logY <- list()
  logYZ1 <- list()
  logYZ2 <- list()

  # Initialize nested pars by batch and source.
  for (n in 1:N){
    Ci [[n]] <- matrix(nrow = P, ncol = K)
    Llist [[n]] <- list()
    Cilist [[n]] <- matrix(NA, nrow = P, ncol = block)
    Lbest [[n]] <- list()
    ProbCi [[n]] <- matrix(0, nrow = P, ncol = K)
    ProbL [[n]] <- list()
    S [[n]] <- list()
    A [[n]] <- list()
    B [[n]] <- list()
    #    Tau[[n]] = list()
    Sigma [[n]] <- list()
    a0 [[n]] <- list()
    b0 [[n]] <- list()
    mu0 [[n]] <- list()
    mu [[n]] <- list()
    L [[n]] <- list()
    Lp [[n]] <- list()
    logC [[n]] <- matrix(nrow = P, ncol = K)
    nu_e [[n]] <- array(1/K, dim = c(P, M[n], K))
    Csize [[n]] <- matrix(nrow = M[n], ncol = K)
    beta [[n]] <- numeric(M[n])
    mpij [[n]] <- matrix(nrow = M[n], ncol = K)

    logY [[n]] <- numeric(M[n])
    logYZ1 [[n]] <- numeric(M[n])
    logYZ2 [[n]] <- numeric(M[n])

    for (m in 1:M[n]){
      Llist [[n]][[m]] <- matrix(NA, nrow = P, ncol = block)
      ProbL [[n]][[m]] <- matrix(0, nrow = P, ncol = K)
      S [[n]][[m]] <- matrix(nrow = D [[n]][[m]], ncol = K)
      A [[n]][[m]] <- matrix(nrow = D [[n]][[m]], ncol = K)
      B [[n]][[m]] <- matrix(nrow = D [[n]][[m]], ncol = K)
      mu [[n]][[m]] <- matrix(nrow = D [[n]][[m]], ncol = K)
      Sigma [[n]][[m]] <- matrix(nrow = D [[n]][[m]], ncol = K)
      StDev <- apply(X [[n]][[m]], 2, sd) #pick a0,b0 based on overall sample variance
      a0 [[n]][[m]] <- rep(1, D [[n]][[m]])
      b0 [[n]][[m]] <- StDev ^ 2
      if (D [[n]][[m]] > 1) mu0 [[n]][[m]] <- colMeans(X [[n]][[m]]) #m0 by overal mean
      if (D [[n]][[m]] == 1) mu0 [[n]][[m]] <- mean(X [[n]][[m]])

      InitL <- cluster::pam(X [[n]][[m]], K) ## Replaced Kmeans for its instability.

      if (m == 1) {now_cluster <- InitL$clustering} else {
        temp_cluster <- clust_align(now_cluster, InitL$clustering, type = 'vec')
        now_cluster <- temp_cluster
      }

      # # initializing cluster centers.
      for (k in 1:K){
        if (sum(InitL$clustering == k) > 1){
          if (D[[n]][[m]] > 1){Sigma [[n]][[m]][, k] <- apply(X [[n]][[m]][InitL$clustering == k, ], MARGIN = 2, FUN = 'sd')
          mu [[n]][[m]][, k] <- apply(X [[n]][[m]][InitL$clustering == k, ], MARGIN = 2, FUN = 'mean')}
          if (D[[n]][[m]] == 1){Sigma [[n]][[m]][, k] <- sd(X [[n]][[m]][InitL$clustering == k, ])
          mu [[n]][[m]][, k] <- mean(X [[n]][[m]][InitL$clustering == k, ])}
        } else {
          if (D [[n]][[m]] > 1){Sigma [[n]][[m]][, k] <- apply(X [[n]][[m]], MARGIN = 2, FUN = 'sd') / K
          mu [[n]][[m]][, k] <- apply(X [[n]][[m]], MARGIN = 2, FUN = 'mean')}
          if (D [[n]][[m]] == 1){Sigma [[n]][[m]][, k] <- sd(X [[n]][[m]] / K)
          mu [[n]][[m]][, k] <- mean(X [[n]][[m]][InitL$clustering == k, ])}
        }
      }

      Lp [[n]][[m]] <- matrix(nrow = P, ncol = K)
      L [[n]][[m]] <- matrix(nrow = P,ncol = K)

    }
  }

  # par in the middle of MCMC updates
  alphaVec <- c() # list of alpha paramters (batch and overall)
  C <- matrix(nrow = P, ncol = K) #overall clustering matrix
  nu_s <- array(1 / K, dim = c(N, P, K))
  Pi <- rep(1 / K, K)
  alpha <- numeric(N)
  if (IndivAlpha) for (n in 1:N) {while (alpha[n] < 1 / K) alpha[n] <- rbeta(1, a, b)}
  if (!IndivAlpha) while (alpha[1] < 1 / K) alpha[] <- rbeta(1, a, b)
  if (IndivBeta) for (n in 1:N){
    for (m in 1:M[n]){while (beta [[n]][m] < 1 / K) beta [[n]][m] <- rbeta(1, a, b)}
  }
  if (!IndivBeta) for(n in 1:N){
    for (m in 1:M[n]){while (beta [[n]][1] < 1 / K) beta [[n]][] = rbeta(1, a, b)}
  }

  record_start <- floor(NumRun/5)
  length_record <- NumRun - record_start

  #Create temp dir to store MCMC outputs
  tmpdir <- paste0('tmp', round(runif(1) * 10000))
  dir.create(tmpdir)
  old_dir <- setwd(tmpdir)
  file.create(paste0('sub', 1:N, '.txt')) #subject clustering files
  ep_files <- paste0(rep(paste0('sub', 1:N), M),
    paste0('ep', rapply(sapply(M, FUN = function(a) c(1:a), simplify = F), identity, how = 'unlist')),'.txt')
  file.create(ep_files)
  file.create(paste0(c('Alpha', 'C'), '.txt'))
  file.create(paste0('Beta', 1:N, '.txt'))

  # Gibbs Sampler
  for (run in 1:NumRun){  # run is the current gibbs iteration
    #Update mixing prob, and ave mixing prob
    for (n in 1:N){
      mpi [n, ] <- margPi(p = Pi, a = alpha[n])
      for (m in 1:M[n]){
        mpij [[n]][m, ] <- margPi(p = mpi [n,],a = beta [[n]][m])
      }
    }

    # Update logF and logL, likelihood given L
    for (n in 1:N){ #batch level
      for (m in 1:M [n]){ #source level
        ll <- matrix(nrow = P, ncol = K)
        logF <- ll
        for (k in 1:K){ # In each clsuter k
          # log-likelihood by cluters
          ll [,k] <- -sum(log(Sigma [[n]][[m]][, k])) + (-D [[n]][[m]] * log(2 * pi) - colSums(((t(X [[n]][[m]]) - mu [[n]][[m]][, k]) / Sigma [[n]][[m]][, k]) ^ 2)) / 2
          logF [,k] <- log(nu_e [[n]][, m, k]) + ll [, k]
        }
        logL <- logF - apply(logF, 1, logSum)
        Lp [[n]][[m]] <- exp(logL)
        ## Updating the expectation for IC calculations.
        if (run >= record_start){
          new.logY <- sum(Lp [[n]][[m]] * ll)
          logY [[n]][[m]] <- (logY [[n]][[m]] * (run - record_start) + new.logY) / (run - record_start + 1)
          new.logYZ2 <- sum(Lp [[n]][[m]] * t(log(mpij [[n]][m,]) + t(ll)))
          logYZ2 [[n]][[m]] <- (logYZ2 [[n]][[m]] * (run - record_start) + new.logYZ2) / (run - record_start + 1)
        }
      }
    }


    #Generate new L and update normal-gamma paramters.
    for (n in 1:N){
      for (m in 1:M [n]){
        for (p in 1:P) L [[n]][[m]][p, ] <- rmultinom(1, 1, Lp[[n]][[m]][p, ])
        if (run > 1) L [[n]][[m]] <- clust_align(C, L [[n]][[m]], type = 'mat') # Avoid label switching by realignment
        Csize [[n]][m, ] <- colSums(L [[n]][[m]])
        ll <- matrix(nrow = P, ncol = K)
        for (k in 1:K){ #NIG updates
          kvec <- X [[n]][[m]][L [[n]][[m]][, k] == 1, ]
          if (D [[n]][[m]] == 1 & Csize [[n]][m, k] > 1){ #1-D source and greater than 1 subjects in k cluster
            S [[n]][[m]][, k] <- sd(kvec) ^ 2
            PostMean <- (sum(kvec) + mu0 [[n]][[m]]) / (Csize [[n]][m, k] + 1)
            B [[n]][[m]][, k] <- b0 [[n]][[m]] + 0.5 * (Csize [[n]][m, k] * S [[n]][[m]][, k] + Csize [[n]][m, k] * (mean(kvec) - mu0 [[n]][[m]]) ^ 2 / (1 + Csize [[n]][m, k]))
          } else if (D [[n]][[m]] > 1 & Csize [[n]][m, k] > 1){ # multi-D source and greater than 1 subjects in k cluster
            PostMean <- (mu0 [[n]][[m]] + rowSums(t(kvec))) /(Csize [[n]][m, k] + 1)
            S [[n]][[m]][, k] <- apply(kvec, MARGIN = 2, FUN = 'sd') ^ 2
            B [[n]][[m]][, k] <- b0 [[n]][[m]] + 0.5 * (Csize [[n]][m, k] * S [[n]][[m]][, k] + Csize [[n]][m, k] * (rowMeans(t(kvec)) - mu0 [[n]][[m]]) ^ 2 / (1 + Csize [[n]][m, k]))
          } else if (Csize [[n]][m, k] == 1){ # Only one subject in cluster k
            PostMean <- (mu0 [[n]][[m]] + t(kvec)) / 2
            B [[n]][[m]][, k] <- b0 [[n]][[m]] + 0.5 * (t(kvec) - mu0 [[n]][[m]]) ^ 2 / 2
          } else {
            PostMean <- mu0 [[n]][[m]]
            B [[n]][[m]][, k] <- b0 [[n]][[m]]
          }
          Lambda <- 1 + Csize [[n]][m, k]
          A [[n]][[m]][, k] <- a0 [[n]][[m]] + Csize [[n]][m, k] / 2
          Tau <- rgamma(D [[n]][[m]], shape = A [[n]][[m]][, k], rate = B [[n]][[m]][, k])
          mu [[n]][[m]][, k] <- rnorm(D [[n]][[m]], PostMean, sd = sqrt(1 / (Tau * Lambda)))
          Sigma [[n]][[m]][, k] <- sqrt(1 / Tau)
          if (A [[n]][[m]][1, k] > 1) {temp.sig <- sqrt(B [[n]][[m]][, k] / (A [[n]][[m]][, k] - 1))}
          else temp.sig <- sqrt(B [[n]][[m]][, k] / A [[n]][[m]][, k])
          ll[, k] <- - sum(log(temp.sig)) + (- D [[n]][[m]] * log(2 * pi) - colSums(((t(X [[n]][[m]])-as.vector(PostMean)) / temp.sig) ^ 2)) / 2
        }
        if (run >= record_start){
          new.logYZ1 <- sum(L [[n]][[m]] * t(as.vector(log(Csize [[n]][m, ] + 1)) - log(P + K) + t(ll)))
          logYZ1 [[n]][[m]] <- (logYZ1[[n]][[m]] * (run - record_start) + new.logYZ1) / (run - record_start + 1)
          ProbL [[n]][[m]] <- (ProbL [[n]][[m]] * (run - record_start) + Lp [[n]][[m]])/(run - record_start + 1)
        }
      }
    }

    #Update betas
    if (run>1){
      for (n in 1:N){
        if (!IndivBeta){
          Numeq <- 0;
          for (m in 1:M[n]) Numeq <- Numeq + sum(L [[n]][[m]] * Ci [[n]])
          for (samp in 1:10){
            betaTemp <- rbeta(1, a + Numeq, b + M[n] * P - Numeq) #generate beta until result>1/K or set to 1/K after 10 tries
            if (betaTemp > 1 / K){
              beta [[n]][] <- betaTemp
              break;
            }
          }
        }
        if (IndivBeta){
          for (m in 1:M[n]){
            Numeq <- sum(L [[n]][[m]] * Ci[[n]])
            for (samp in 1:10){
              betaTemp <- rbeta(1, a + Numeq, b + P - Numeq)
              if (betaTemp > 1 / K){
                beta [[n]][m] <- betaTemp
                break;
              }
            }
          }
        }
      }
    }

    #Update Ci's
    for (n in 1:N){
      logC [[n]] <- log(nu_s[n, , ])
      for (m in 1:M[n]){
        logC[[n]] <- logC[[n]] + L[[n]][[m]] * log(beta[[n]][m]) + (1-L[[n]][[m]]) * (log(1 - beta[[n]][m]) - log(K - 1))
      }
      logC[[n]] <- logC[[n]] - apply(logC[[n]], 1, logSum)
    }
    LC <- lapply(logC,exp)
    for (n in 1:N){
      for (p in 1:P){
        Ci[[n]][p, ] <- rmultinom(1, 1, (LC [[n]][p, ] + adj))
      }
      if (run > 1) Ci[[n]] <- clust_align(C, Ci[[n]], type = 'mat')
      if (run >= record_start) ProbCi[[n]] <- (ProbCi[[n]] * (run-record_start) + LC[[n]])/(run - record_start + 1)
    }
    #Update alphas
    if (run > 1){
      if (!IndivAlpha){
        Numeq <- 0;
        for (n in 1:N) Numeq <- Numeq + sum(Ci[[n]] * C)
        for (count in 1:10) {AlphaTemp <- rbeta(1, a + Numeq, b + N * P - Numeq)
        if (AlphaTemp > 1 / K){
          alpha[] <- AlphaTemp
          break;}}}
      if (IndivAlpha){
        for (n in 1:N){
          Numeq <- sum(Ci[[n]] * C)
          for (count in 1:10) {AlphaTemp <- rbeta(1, a + Numeq, b + P - Numeq)
          if (AlphaTemp > 1 / K){
            alpha[n] <- AlphaTemp
            break;}}
        }
      }
    }
    #Update C's
    Cprob <- matrix(log(Pi), nrow = P, ncol = K, byrow = T)
    for (n in 1:N){
      Cprob <- Cprob + Ci[[n]] * log(alpha[n]) + (1-Ci[[n]]) * (log(1 - alpha[n]) - log(K - 1))
    }
    Cprob <- Cprob - apply(Cprob, 1, logSum)
    Cprob <- exp(Cprob)
    for(p in 1:P){
      C[p, ] <- rmultinom(1, 1, Cprob[p, ])
    }
    if(run >= record_start) ProbC <- (ProbC * (run-record_start) + Cprob)/(run-record_start + 1)

    #Update Pi
    Gamma <- Concentration+colSums(C)
    Pi <- gtools::rdirichlet(1, Gamma)

    update_index <- (run - 1) %% block + 1

    if (IndivAlpha) alphaVec <- cbind(alphaVec, alpha)
    if (!IndivAlpha) alphaVec[update_index] <- alpha[1]
    if (IndivBeta){for (n in 1:N){Betalist[[n]] <- cbind(Betalist[[n]], beta[[n]])}}
    if (!IndivBeta){for (n in 1:N){Betalist[[n]][update_index] <- beta[[n]][1]}}

    #Update nu_e and nu_s
    for (n in 1:N){
      for (m in 1:M[n]){
        nu_e [[n]][, m, ] <- Ci[[n]] * beta[[n]][m] + (1 - Ci[[n]]) * (1 - beta[[n]][m]) / (K - 1)
      }
      nu_s[n, , ] <- C * alpha[n] + (1 - C) * (1 - alpha[n]) / (K - 1)
    }

    #update lists of L,Ci,C estimates of current Run
    Clist[, update_index] <- C %*% c(1:K)
    for (n in 1:N){
      Cilist[[n]][, update_index] <- Ci[[n]] %*% c(1:K)
      for (m in 1:M[n]){
        Llist[[n]][[m]][, update_index] <- L[[n]][[m]] %*% c(1:K)
        if(run %% block == 0| run == NumRun){
          filename <- paste0('sub',n,'ep',m,'.txt')
          write(t(na.omit(t(Llist[[n]][[m]]))), file = filename, ncolumns = P, append = T, sep = ' ')
          Llist[[n]][[m]] <- matrix(NA,nrow = P,ncol = block)
        }
      }
      if (run %% block == 0 | run == NumRun){
        filename <- paste0('sub',n,'.txt')
        write(t(na.omit(t(Cilist[[n]]))), file = filename, ncolumns = P, append = T, sep = ' ')
        Cilist[[n]] <- matrix(NA,nrow = P, ncol = block)
        filename <- paste0('Beta',n,'.txt')
        if (IndivBeta){
          write(na.omit(Betalist[[n]]), file = filename, ncolumns = M[n], append = T, sep = ' ')
        } else{
          write(na.omit(Betalist[[n]]), file = filename, ncolumns = 1, append = T, sep = ' ')
        }
      }
    }
    if (run %% block ==0| run == NumRun){
      write(t(na.omit(t(Clist))), file = 'C.txt', ncolumns = P, append = T, sep = ' ')
      Clist <- matrix(NA, nrow = P, ncol = block)
      if (IndivAlpha){
        write(na.omit(alphaVec), file = 'Alpha.txt', ncolumns = N, append = T, sep = ' ')
      } else{
        write(na.omit(alphaVec), file = 'Alpha.txt', ncolumns = 1, append = T, sep = ' ')
      }
      alphaVec <- c()
      cat('Elapsed',difftime(Sys.time(),time.start, units = 'hours'),'hr,',
        difftime(Sys.time(),time.start, units = 'hours')/run*(NumRun-run),'hr to go','\n')
      Betalist <- vector("list", N)
    }
  }

  # Estimation from posterior samples.
  # # Output Alphas and Betas.
  alpha <- read.table('Alpha.txt', sep = ' ', skip = floor(NumRun / 5))
  colnames(alpha) <- paste0('alpha', 1:dim(alpha)[2])
  Alpha <- apply(alpha,MARGIN = 2, quantile, probs = 0.5)
  AlphaBounds <- cbind(apply(alpha,MARGIN = 2, quantile, probs = 0.025),
    apply(alpha,MARGIN = 2, quantile, probs = 0.975))
  colnames(AlphaBounds) <- c('lower','upper')
  rm(alpha)

  BetaBounds <- list()
  Betas <- list()
  for (n in 1:N){
    beta <- read.table(paste0('Beta',n,'.txt'), sep = ' ', skip = floor(NumRun / 5))
    colnames(beta) <- paste0('beta',1:dim(beta)[2])
    Betas[[n]] <- apply(beta,MARGIN = 2, quantile, probs = 0.5)
    BetaBounds[[n]] <- cbind(apply(beta,MARGIN = 2, quantile, probs = 0.025),
      apply(beta,MARGIN = 2, quantile, probs = 0.975))
    colnames(BetaBounds[[n]]) <- c('lower','upper')
    rm(beta)
  }
  # # Choose hard clustering(point estimate) by least square as in (Dahl,2006)
  labellist <- as.matrix(read.table('C.txt',sep = ' ', skip = floor(NumRun / 5)))
  Cout <- HardCluster(labellist)
  Cbest <- Cout$Cbest
  Cvar <- mean(Cout$vec)
  Cvec <- Cout$vec
  Civar <- numeric(N)
  Civec <- list()
  Lvar <- list()

  for (n in 1:N){
    labellist <- as.matrix(read.table(paste0('sub',n,'.txt'), sep = ' ', skip = floor(NumRun / 5)))
    Cout <- HardCluster(labellist)
    Cibest[,n] <- Cout$Cbest
    Civar[n] <- mean(Cout$vec)
    Civec[[n]] <- Cout$vec
    Lvar[[n]] <- numeric(M[n])

    for (m in 1:M[n]){
      labellist <- as.matrix(read.table(paste0('sub',n,'ep',m,'.txt'), sep = ' ', skip = floor(NumRun / 5)))
      Cout <- HardCluster(labellist)
      Lbest[[n]][[m]] <- Cout$Cbest
      Lvar[[n]][m] <- mean(Cout$vec)
    }
  }

  # # Information Output
  llY <- list()
  L <- rapply(Lbest, function(x) t(sapply(x, function(a) a == c(1:K))+0), how = 'list')
  for (n in 1:N){
    llY[[n]] <- list()
    for (m in 1:M[n]){
      Csize[[n]][m,] <- colSums(L[[n]][[m]])
      ll <- matrix(nrow = P, ncol = K)
      for (k in 1:K){ #NIG updates
        kvec <- X[[n]][[m]][L[[n]][[m]][, k] == 1, ]
        if (D[[n]][[m]] == 1 & Csize[[n]][m, k] > 1){ #1-D source and greater than 1 subjects in k cluster
          S[[n]][[m]][, k] <- sd(kvec)^2
          mu[[n]][[m]][, k] <- (sum(kvec) + mu0[[n]][[m]]) / (Csize[[n]][m, k] + 1)
          B[[n]][[m]][,k] <- b0[[n]][[m]] + 0.5 * (Csize[[n]][m, k] * S[[n]][[m]][, k] + Csize[[n]][m, k] * (mean(kvec) - mu0[[n]][[m]]) ^ 2 / (1 + Csize[[n]][m, k]))
        } else if (D[[n]][[m]] > 1 & Csize[[n]][m, k] > 1){ # multi-D source and greater than 1 subjects in k cluster
          mu[[n]][[m]][, k] <- (mu0[[n]][[m]] + rowSums(t(kvec))) / (Csize[[n]][m, k] + 1)
          S[[n]][[m]][, k] <- apply(kvec, MARGIN = 2, FUN = 'sd') ^ 2
          B[[n]][[m]][, k] <- b0[[n]][[m]] + 0.5 * (Csize[[n]][m, k] * S[[n]][[m]][, k] + Csize[[n]][m, k] * (rowMeans(t(kvec)) - mu0[[n]][[m]]) ^ 2 / (1 + Csize[[n]][m, k]))
        } else if (Csize[[n]][m, k] == 1){ # Only one subject in cluster k
          mu[[n]][[m]][, k] <- (mu0[[n]][[m]] + t(kvec)) / 2
          B[[n]][[m]][, k] <- b0[[n]][[m]] + 0.5 * (t(kvec) - mu0[[n]][[m]]) ^ 2 / 2
        } else {
          mu[[n]][[m]][, k] <- mu0[[n]][[m]]
          B[[n]][[m]][, k] <- b0[[n]][[m]]
        }
        A[[n]][[m]][, k] <- a0[[n]][[m]] + Csize[[n]][m, k] / 2
        if (A[[n]][[m]][1, k] > 1) {Sigma[[n]][[m]][, k] <- sqrt(B[[n]][[m]][, k] / (A[[n]][[m]][, k] - 1))}
        else Sigma[[n]][[m]][, k] <- sqrt(B[[n]][[m]][, k] / A[[n]][[m]][, k])

        ll[, k] <- -sum(log(Sigma[[n]][[m]][, k])) + ( -D[[n]][[m]] * log(2 * pi) - colSums((((t(X[[n]][[m]]) - mu[[n]][[m]][, k]) / Sigma[[n]][[m]][, k]) ^ 2))) / 2
      }
      llY[[n]][[m]] <- sum(ll * L[[n]][[m]])
    }
  }

  ## Removign the temp files and dir
  setwd(old_dir)
  unlink(tmpdir, recursive = T)

  BIC <- 2 * sum(rapply(llY, sum, how = 'unlist')) / Total_e - (2 * K * sum(rapply(D, sum, how = 'unlist')) + K - 1) * log(P) / Total_e
  BIC2 <- 2 * sum(rapply(logYZ1, sum, how = 'unlist')) / Total_e - (2 * K * sum(rapply(D, sum, how = 'unlist')) + K - 1) * log(P) / Total_e
  DIC4 <- - sum(rapply(logYZ1, sum, how = 'unlist')) / Total_e + 2 * sum(rapply(logYZ2, sum, how = 'unlist')) / Total_e
  DIC7 <- - sum(rapply(llY, sum, how = 'unlist'))/Total_e + 2 * sum(rapply(logY, sum, how = 'unlist')) / Total_e
  COH <- (K * mean(Alpha) - 1) / (K - 1)
  # Output Organization
  return(list(Alpha = Alpha,
              AlphaBounds = AlphaBounds,
              Betas = Betas,
              BetaBounds = BetaBounds,
              Cbest = Cbest,
              Cibest = Cibest,
              Lbest = Lbest,
              ICs = list(BIC,BIC2,DIC4,DIC7,COH),
              Cvar = Cvar,
              Civar = Civar,
              Lvar = Lvar,
              Cvec = Cvec,
              Civec = Civec,
              ProbC = ProbC,
              ProbCi = ProbCi,
              ProbL = ProbL))
}