"bsad" <- function(x, xmin, xmax, nint, MaxNCos, mcmc = list(), prior = list(), smoother = c("geometric", "algebraic"), 
                   parametric = c("none", "normal", "gamma", "laplace"), marginal.likelihood = TRUE) {
  cl <- match.call()
  xdata <- x
  nobs <- length(xdata)
  if (missing(xmin)) {
    xmin <- min(xdata) - sd(xdata)/2
  }
  if (missing(xmax)) {
    xmax <- max(xdata) + sd(xdata)/2
  }
  if (missing(nint)) {
    nint <- 201
  }

  xdelta <- (xmax - xmin)/nint
  xgrid <- seq(xmin + xdelta/2, by = xdelta, length = nint)
  xgridm <- mean(xgrid)

  xgl <- xgrid - xdelta/2
  xgu <- xgrid + xdelta/2
  xdatamat <- matrix(xdata, nrow = nobs, ncol = nint)
  ndata <- apply(sweep(xdatamat, 2, xgl, FUN = ">") * sweep(xdatamat, 2, xgu, FUN = "<="), 2, sum)

  if (missing(MaxNCos)) {
    MaxNCos <- min(nobs, nint) - 2
  } else {
    nbasis <- min(nobs, nint) - 2
    if (nbasis < MaxNCos) {
      stop("The maximun number of Fourier coefficents should be less than
           \t\t\t      the minimum number of observations and gridpoints minus 2.")
    }
    }

  smoother <- match.arg(smoother)
  if (smoother == "geometric")
    smooth <- 1 else smooth <- 0

  if (smooth == 1) {
    kall <- 1:MaxNCos
  } else {
    kall <- log(1:MaxNCos)
  }

  if (MaxNCos >= 20) {
    probk <- matrix(0, nrow = 19, ncol = 19)
    for (i in 1:nrow(probk)) {
      for (j in 1:ncol(probk)) {
        if (i <= j) {
          probk[i, j] <- ((10 + i) - j) * ((10 + i) > j)
          probk[j, i] <- probk[i, j]
        }
      }
    }
  } else if (MaxNCos < 20 && MaxNCos >= 8) {
    probk <- matrix(0, nrow = 7, ncol = 7)
    for (i in 1:nrow(probk)) {
      for (j in 1:ncol(probk)) {
        if (i <= j) {
          probk[i, j] <- ((4 + i) - j) * ((4 + i) > j)
          probk[j, i] <- probk[i, j]
        }
      }
    }
  } else {
    probk <- matrix(0, nrow = 3, ncol = 3)
    probk[1, ] <- c(3, 1, 0)
    probk[2, ] <- c(1, 3, 1)
    probk[3, ] <- c(1, 3, 0)
  }

  probk <- sweep(probk, 2, colSums(probk), FUN = "/")
  cdfk <- apply(probk, 2, cumsum)
  Probzz <- 0.5 * matrix(1, nrow = 2, ncol = 1)

  parametric <- match.arg(parametric)
  if (parametric == "none") {
    flagnormal <- 0
    flaggamma <- 0
    flaglaplace <- 0
  }
  if (parametric == "normal") {
    flagnormal <- 1
    flaggamma <- 0
    flaglaplace <- 0
  }
  if (parametric == "gamma") {
    flagnormal <- 0
    flaggamma <- 1
    flaglaplace <- 0
  }
  if (parametric == "laplace") {
    flagnormal <- 0
    flaggamma <- 0
    flaglaplace <- 1
  }

  d <- intsim(rep(1, nint), xdelta)
  if (flagnormal == 1) {
    cx <- intsim((xgrid - xgridm)^2, xdelta)
    cx2 <- cx/d
    dmat <- cbind((xgrid - xgridm), ((xgrid - xgridm)^2 - cx2))
  }

  if (flaggamma == 1) {
    dmat <- log(xgrid - xmin)
    clnx <- mean(dmat)
    dmat <- dmat - clnx
    dmat <- cbind(dmat, (xgrid - xgridm))
  }
  if (flaglaplace == 1) {
    cx <- intsim(abs(xgrid - xgridm), xdelta)
    cxabs <- cx/d
    dmat <- cbind(abs(xgrid - xgridm) - cxabs)
  }

  if (flagnormal == 0 && flaggamma == 0 && flaglaplace == 0) {
    dmat <- matrix(1, nrow = nint, ncol = 1)
  }
  dtd <- crossprod(dmat)
  npar <- ncol(dmat)

  ugrid <- (2 * (1:nint) - 1)/(2 * nint)

  phi <- sqrt(2) * cos(pi * tcrossprod(ugrid, kall))

  phidata <- sqrt(2) * cos(pi * tcrossprod(((xdata - xmin)/(xmax - xmin)), kall))

  if (flagnormal == 1) {
    ddata <- cbind((xdata - xgridm), ((xdata - xgridm)^2 - cx2))
  }

  if (flaggamma == 1) {
    ddata <- log(xdata - xmin) - clnx
    ddata <- cbind(ddata, (xdata - xgridm))
  }

  if (flaglaplace == 1) {
    ddata <- cbind(abs(xdata - xgridm) - cxabs)
  }

  if (flagnormal == 0 && flaggamma == 0 && flaglaplace == 0) {
    ddata <- matrix(1, nrow = nobs, ncol = 1)
  }

  privals <- list(gmax = 5, PriorProbs = c(0.5, 0.5), beta_m0 = matrix(0, nrow = npar, ncol = 1), beta_v0 = 10 * diag(npar),
                  r0 = 4, s0 = 4, u0 = 4, v0 = 4, PriorKappa = matrix(1/(MaxNCos + 1), nrow = MaxNCos + 1, ncol = 1), KappaGrid = 0:MaxNCos)
  privals[names(prior)] <- prior

  gmax <- privals$gmax

  PriorProbs <- privals$PriorProbs
  nmodels <- length(PriorProbs)
  modnames <- c("Para", "SemiPar")

  BetaMean0 <- privals$beta_m0
  BetaVar0 <- privals$beta_v0
  BetaVari0 <- solve(BetaVar0)
  BetaVariMean0 <- BetaVari0 %*% BetaMean0

  r0 <- privals$r0
  s0 <- privals$s0

  u0 <- privals$u0
  v0 <- privals$v0

  PriorKappa <- privals$PriorKappa
  KappaGrid <- privals$KappaGrid
  cdfPriorKappa <- apply(PriorKappa, 2, cumsum)
  lnPriorKappa <- log(PriorKappa)

  mcvals <- list(nblow = 10000, smcmc = 1000, nskip = 10, ndisp = 1000, kappaloop = 5)
  mcvals[names(mcmc)] <- mcmc
  nblow <- mcvals$nblow
  smcmc <- mcvals$smcmc
  nskip <- mcvals$nskip
  ndisp <- mcvals$ndisp
  kappaloop <- mcvals$kappaloop

  mcmctime <- system.time({
    foo <- .Fortran("bsad", as.integer(nblow), as.integer(smcmc), as.integer(nskip),
                    as.integer(ndisp), as.integer(kappaloop),
                    as.integer(nobs), as.integer(nint), as.integer(npar), as.integer(MaxNCos),
                    as.matrix(probk), as.matrix(cdfk), as.integer(ncol(probk)),
                    as.double(gmax), as.integer(smooth), as.double(PriorProbs),
                    as.integer(nmodels), as.matrix(BetaVari0), as.double(BetaVariMean0),
                    as.double(r0), as.double(s0), as.double(u0), as.double(v0),
                    as.matrix(dmat), as.matrix(ddata), as.matrix(dtd),
                    as.matrix(phi), as.matrix(phidata), as.integer(ndata),
                    as.double(xdelta), as.double(lnPriorKappa), as.double(cdfPriorKappa),
                    as.integer(KappaGrid), betaParg = matrix(0, nrow = smcmc, ncol = npar),
                    sigmaParg = numeric(smcmc), betag = matrix(0, nrow = smcmc, ncol = npar),
                    sigmag = numeric(smcmc), taug = numeric(smcmc), gamg = numeric(smcmc),
                    thetag = matrix(0,  nrow = smcmc, ncol = MaxNCos), kappag = as.integer(numeric(smcmc)),
                    betaMaxKappag = matrix(0, nrow = smcmc, ncol = npar),
                    sigmaMaxKappag = numeric(smcmc), tauMaxKappag = numeric(smcmc),
                    gamMaxKappag = numeric(smcmc), thetaMaxKappag = matrix(0, nrow = smcmc, ncol = MaxNCos),
                    fparg = matrix(0, nrow = smcmc, ncol = nint), fsemig = matrix(0, nrow = smcmc, ncol = nint),
                    fsemiMaxKappag = matrix(0, nrow = smcmc, ncol = nint), Probzzg = matrix(0, nrow = smcmc, ncol = nmodels),
                    ExactLogLikeg = matrix(0, nrow = smcmc, ncol = 3), MetProbPar = as.integer(numeric(1)),
                    MetProbSemi = as.integer(numeric(1)), NAOK = TRUE, PACKAGE = "bsamGP")
  })

  mcmc.draws <- list()
  mcmc.draws$betaPar <- foo$betaParg
  mcmc.draws$sigmaPar <- foo$sigmaParg
  mcmc.draws$beta <- foo$betag
  mcmc.draws$sigma <- foo$sigmag
  mcmc.draws$tau <- foo$taug
  mcmc.draws$gam <- foo$gamg
  mcmc.draws$theta <- foo$thetag
  mcmc.draws$kappa <- foo$kappag
  mcmc.draws$betaMaxKappa <- foo$betaMaxKappag
  mcmc.draws$sigmaMaxKappa <- foo$sigmaMaxKappag
  mcmc.draws$tauMaxKappa <- foo$tauMaxKappag
  mcmc.draws$gamMaxKappa <- foo$gamMaxKappag
  mcmc.draws$thetaMaxKappa <- foo$thetaMaxKappag
  mcmc.draws$Probzz <- foo$Probzzg

  fit.draws <- list()
  fit.draws$fpar <- foo$fparg
  fit.draws$fsemi <- foo$fsemig
  fit.draws$fsemiMaxKappa <- foo$fsemiMaxKappag

  post.est <- list()
  post.est$betaParm <- apply(foo$betaParg, 2, mean)
  post.est$betaPars <- apply(foo$betaParg, 2, sd)
  post.est$sigmaParm <- mean(foo$sigmaParg)
  post.est$sigmaPars <- sd(foo$sigmaParg)

  post.est$betam <- apply(foo$betag, 2, mean)
  post.est$betas <- apply(foo$betag, 2, sd)
  post.est$sigmam <- mean(foo$sigmag)
  post.est$sigmas <- sd(foo$sigmag)
  post.est$taum <- mean(foo$taug)
  post.est$taus <- sd(foo$taug)
  post.est$gamm <- mean(foo$gamg)
  post.est$gams <- sd(foo$gamg)
  post.est$kappam <- mean(foo$kappag)
  post.est$kappas <- sd(foo$kappag)
  thetam <- apply(foo$thetag, 2, mean)
  post.est$thetam <- thetam
  thetas <- apply(foo$thetag, 2, function(x) x^2)
  post.est$thetas <- sqrt(abs(thetas - smcmc * thetam^2)/smcmc)

  post.est$betaMaxKappam <- apply(foo$betaMaxKappag, 2, mean)
  post.est$betaMaxKappas <- apply(foo$betaMaxKappag, 2, sd)
  post.est$sigmaMaxKappam <- mean(foo$sigmaMaxKappag)
  post.est$sigmaMaxKappas <- sd(foo$sigmaMaxKappag)
  post.est$tauMaxKappam <- mean(foo$tauMaxKappag)
  post.est$tauMaxKappas <- sd(foo$tauMaxKappag)
  post.est$gamMaxKappam <- mean(foo$gamMaxKappag)
  post.est$gamMaxKappas <- sd(foo$gamMaxKappag)
  thetaMaxKappam <- apply(foo$thetaMaxKappag, 2, mean)
  post.est$thetaMaxKappam <- thetaMaxKappam
  thetaMaxKappas <- apply(foo$thetaMaxKappag, 2, function(x) x^2)
  post.est$thetaMaxKappas <- thetaMaxKappas
  post.est$thetaMaxKappas <- sqrt(abs(thetaMaxKappas - smcmc * thetaMaxKappam^2)/smcmc)

  out <- list()
  out$call <- cl
  out$mcmctime <- mcmctime

  out$x <- xdata
  out$xdelta <- xdelta
  out$xgrid <- xgrid
  out$dmat <- dmat
  out$ddata <- ddata
  out$nobs <- nobs
  out$nint <- nint
  out$xmin <- xmin
  out$xmax <- xmax
  out$MaxNCos <- MaxNCos

  out$mcmc <- mcvals
  out$prior <- privals

  out$parametric <- parametric
  out$smooth <- smoother

  out$mcmc.draws <- mcmc.draws
  out$fit.draws <- fit.draws
  out$loglike.draws <- foo$ExactLogLikeg

  out$post.est <- post.est

  out$MetProb <- c(foo$MetProbPar/smcmc, foo$MetProbSemi/smcmc)
  names(out$MetProb) <- c("Param", "Semiparam")

  Probzzm <- apply(foo$Probzzg, 2, mean)
  out$PostProbs <- (1 - rev(Probzzm))/(2 - sum(Probzzm))

  out$marginal.likelihood <- marginal.likelihood
  if (marginal.likelihood) {
    out$maxloglike <- apply(foo$ExactLogLikeg, 2, max)
    out$lmarg <- out$maxloglike - log(apply(exp(-sweep(foo$ExactLogLikeg, 2, out$maxloglike)), 2, mean))
    names(out$lmarg) <- c("Param", "Semiparam", "SemiparamMaxKappa")

    if (flagnormal == 0 && flaggamma == 0 && flaglaplace == 0) {
      out$maxloglike <- out$maxloglike[-1]
      out$lmarg <- out$lmarg[-1]
      names(out$lmarg) <- c("Semiparam", "SemiparamMaxKappa")
    }
  }

  class(out) <- "bsad"
  out
}
