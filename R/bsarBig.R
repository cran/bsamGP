"bsarBig" <- function(formula, nbasis, nint, mcmc = list(), prior = list(), verbose = FALSE) {
  cl <- match.call()

  yxdata <- parse.formula(formula, data = NULL)
  ydata <- yxdata[[1]]
  yname <- yxdata[[2]]
  xdata <- yxdata[[3]][, -1, drop = FALSE]
  xname <- yxdata[[4]][-1]
  nobs <- length(ydata)

  if (missing(nbasis)) {
    stop("The number of basis functions is specified by user.")
  }

  if (missing(nint)) {
    nint <- 500
  }

  xdelta <- 1/nint
  xgrid <- seq(xdelta/2, 1 - xdelta/2, xdelta)

  privals <- list(tau2_m0 = 1, tau2_v0 = 100, w0 = 2, sigma2_r0 = 5, sigma2_s0 = 5)
  privals[names(prior)] <- prior
  w0 <- privals$w0
  tau2_m0 <- privals$tau2_m0
  tau2_v0 <- privals$tau2_v0
  sigma2_r0 <- privals$sigma2_r0
  sigma2_s0 <- privals$sigma2_s0

  mcvals <- list(nblow = 1000, smcmc = 1000, nskip = 1, ndisp = 500)
  mcvals[names(mcmc)] <- mcmc
  nblow <- mcvals$nblow
  smcmc <- mcvals$smcmc
  nskip <- mcvals$nskip
  ndisp <- mcvals$ndisp

  mcmctime <- system.time({
    foo <- .Fortran("bsarbig", as.integer(verbose), as.double(ydata), as.double(xdata),
                    as.integer(nobs), as.integer(nint), as.integer(nbasis), as.double(w0),
                    as.double(tau2_m0), as.double(tau2_v0), as.double(sigma2_r0), as.double(sigma2_s0),
                    as.integer(nblow), as.integer(smcmc), as.integer(nskip), as.integer(ndisp),
                    fhatg = matrix(0, nrow = smcmc, ncol = nint), thetag = matrix(0, nrow = smcmc, ncol = nbasis),
                    sigmag = numeric(smcmc), smoothg = matrix(0, nrow = smcmc, ncol = 2),
                    NAOK = TRUE, PACKAGE = "bsamGP")
  })
  mcmc.draws <- list()
  mcmc.draws$tau2 <- foo$smoothg[,1]^2
  mcmc.draws$tau <- foo$smoothg[,1]
  mcmc.draws$gamma <- foo$smoothg[,2]
  mcmc.draws$theta <- foo$thetag
  mcmc.draws$sigma2 <- foo$sigmag^2
  mcmc.draws$sigma <- foo$sigmag

  fit.draws <- list()
  fit.draws$xgrid <- xgrid
  fit.draws$fxgrid <- foo$fhatg

  # Posterior Mean
  fhatm <- apply(foo$fhatg, 2, mean)
  thetam <- apply(foo$thetag, 2, mean)
  sigmam <- mean(foo$sigmag)
  taum <- mean(foo$smoothg[,1])
  gammam <- mean(foo$smoothg[,2])

  # Posterior STD DEV
  fhats <- apply(foo$fhatg, 2, sd)
  thetas <- apply(foo$thetag, 2, sd)
  sigmas <- sd(foo$sigmag)
  taus <- sd(foo$smoothg[,1])
  gammas <- sd(foo$smoothg[,2])

  # Posterior Quantiles
  fhatq <- t(apply(foo$fhatg, 2, quantile, probs=c(.5,.025,.975)))
  thetaq <- t(apply(foo$thetag, 2, quantile, probs=c(.5,.025,.975)))
  sigmaq <- quantile(foo$sigmag, probs=c(.5,.025,.975))
  tauq <- quantile(foo$smoothg[,1], probs=c(.5,0.25,.975))
  gammaq <- quantile(foo$smoothg[,2], probs=c(.5,0.25,.975))

  post.est <- list()
  post.est$fhatm <- fhatm
  post.est$fhats <- fhats
  post.est$fhatq <- fhatq
  post.est$thetam <- thetam
  post.est$thetas <- thetas
  post.est$thetaq <- thetaq
  post.est$sigmam <- sigmam
  post.est$sigmas <- sigmas
  post.est$sigmaq <- sigmaq
  post.est$taum <- taum
  post.est$taus <- taus
  post.est$tauq <- tauq
  post.est$gammam <- gammam
  post.est$gammas <- gammas
  post.est$gammaq <- gammaq

  res.out <- list()
  res.out$call <- cl
  res.out$model <- 'bsarBig'
  res.out$y <- ydata
  res.out$x <- xdata
  res.out$n <- nobs
  res.out$nint <- nint
  res.out$nbasis <- nbasis

  res.out$prior <- privals

  res.out$mcmctime <- mcmctime
  res.out$mcmc <- mcvals

  res.out$mcmc.draws <- mcmc.draws
  res.out$fit.draws <- fit.draws

  res.out$post.est <- post.est

  class(res.out) <- "bsam"
  res.out
}
