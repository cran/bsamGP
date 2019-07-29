"bsaq" <- function(formula, xmin, xmax, p, nbasis, nint, mcmc = list(), prior = list(),
                   shape = c("Free", "Increasing", "Decreasing", "IncreasingConvex",
                             "DecreasingConcave", "IncreasingConcave", "DecreasingConvex",
                             "IncreasingS", "DecreasingS", "IncreasingRotatedS",
                             "DecreasingRotatedS", "InvertedU", "Ushape","IncMultExtreme","DecMultExtreme"), nExtreme = NULL,
                   marginal.likelihood = TRUE, spm.adequacy = FALSE, verbose = FALSE) {
  cl <- match.call()

  if (missing(p)) {
    p <- 0.5
  } else {
    if (p <= 0 || p >= 1) {
      stop("p must be in (0,1).\n")
    }
  }

  ywxdata <- interpret.bsam(formula)
  yobs <- ywxdata[[1]]
  yname <- ywxdata[[2]]
  wdata <- ywxdata[[3]]
  wnames <- ywxdata[[4]]
  xobs <- ywxdata[[5]]
  xname <- ywxdata[[6]]

  nobs <- nrow(yobs)
  nparw <- ncol(wdata)
  ndimw <- nparw - 1

  if (missing(nbasis))
    stop("The number of basis functions are specified by user.")

  fshape <- function.shape(shape)
  fmodel <- fshape$fmodel
  fpm <- fshape$fpm
  nfun <- fshape$nfun

  if (nfun != ncol(xobs))
    stop("The number of shape and columns of x should be same.")

  if (missing(nint)) {
    nint <- 200
  }
  if (missing(xmin)) {
    xmin <- apply(xobs, 2, min)
  } else {
    if (nfun != length(xmin))
      stop("The number of shape and length of xmin should be same.")
  }
  if (missing(xmax)) {
    xmax <- apply(xobs, 2, max)
  } else {
    if (nfun != length(xmax))
      stop("The number of shape and length of xmax should be same.")
  }
  xmid <- (xmin + xmax)/2
  xrange <- xmax - xmin
  xdelta <- (xrange)/nint
  xgrid <- matrix(0, nrow = nint + 1, ncol = nfun)
  for (i in 1:nfun) {
    xgrid[, i] <- seq(xmin[i], by = xdelta[i], length = nint + 1)
  }

  if(length(nExtreme) != length(which(fmodel==8)))
    stop("The length of extreme points are less than specified multiple extreme shape constraints.")
  nExtremes <- numeric(nfun)
  nExtremes[fmodel==7] <- 1
  nExtremes[fmodel==8] <- nExtreme

  if (any((fmodel==8) & (nExtremes==1))) {
      warning("MultiExtrema with 1 extreme point is equal to either U shape or inverted U shape.")
      fmodel[(fmodel==8) & (nExtremes==1)] <- 7
  }

  maxNext = ifelse(max(nExtremes)>1, max(nExtremes), 1)
  if(maxNext > 1) {
    mPos = which(nExtremes!=0)
    omega_m0 <- matrix(NA, nrow=maxNext, ncol=nfun)
    for(pos in mPos)
      omega_m0[1:nExtremes[pos], pos] <- seq(xmin[pos], xmax[pos], length.out=nExtremes[pos]+2)[-c(1,nExtremes[pos]+2)]
  } else {
    omega_m0 <- matrix((xmin + xmax)/2, nrow=1, ncol=nfun)
  }

  privals <- list(iflagprior = 0, theta0_m0 = 0, theta0_s0 = 100, tau2_m0 = 1, tau2_v0 = 100, w0 = 2, beta_m0 = numeric(nparw),
                  beta_v0 = diag(100, nparw), sigma2_m0 = 1, sigma2_v0 = 1000, alpha_m0 = 3, alpha_s0 = 50, iflagpsi = 1, psifixed = 100,
                  omega_m0 = omega_m0, omega_s0 = xrange/8)
  privals[names(prior)] <- prior

  iflagprior <- privals$iflagprior
  theta0_m0 <- privals$theta0_m0
  theta0_s0 <- privals$theta0_s0
  tau2_m0 <- privals$tau2_m0
  tau2_v0 <- privals$tau2_v0
  w0 <- privals$w0
  beta_m0 <- privals$beta_m0
  beta_v0 <- privals$beta_v0
  sigma2_m0 <- privals$sigma2_m0
  sigma2_v0 <- privals$sigma2_v0
  alpha_m0 <- privals$alpha_m0
  alpha_s0 <- privals$alpha_s0
  iflagpsi <- privals$iflagpsi
  psifixed <- privals$psifixed
  psi_m0 <- psifixed
  psi_s0 <- 10 * psifixed
  omega_m0 <- privals$omega_m0
  omega_s0 <- privals$omega_s0

  mcvals <- list(nblow0 = 1000, nblow = 10000, smcmc = 1000, nskip = 10, ndisp = 1000, maxmodmet = 5)
  mcvals[names(mcmc)] <- mcmc
  nblow0 <- mcvals$nblow0
  nblow <- mcvals$nblow
  smcmc <- mcvals$smcmc
  nskip <- mcvals$nskip
  ndisp <- mcvals$ndisp
  maxmodmet <- mcvals$maxmodmet

  if (max(fmodel) == 1)
    maxmodmet <- 0

  mcmctime <- system.time({
    foo <- .Fortran("bsaqam", as.integer(verbose), as.double(yobs), as.matrix(wdata), as.matrix(xobs), as.integer(nobs), as.integer(nparw),
                    as.integer(nfun), as.integer(nbasis), as.integer(nint), as.integer(nExtremes), as.integer(maxNext), as.integer(fmodel), as.double(fpm), as.double(p), as.double(theta0_m0),
                    as.double(theta0_s0), as.double(tau2_m0), as.double(tau2_v0), as.double(w0), as.double(beta_m0), as.matrix(beta_v0),
                    as.double(alpha_m0), as.double(alpha_s0), as.double(psi_m0), as.double(psi_s0), as.double(psifixed), as.matrix(omega_m0),
                    as.double(omega_s0), as.double(sigma2_m0), as.double(sigma2_v0), as.integer(iflagprior), as.integer(iflagpsi),
                    as.integer(maxmodmet), as.integer(nblow0), as.integer(nblow), as.integer(smcmc), as.integer(nskip), as.integer(ndisp),
                    zetag = matrix(0, smcmc, nfun), tau2g = matrix(0, smcmc, nfun), gammag = matrix(0, smcmc, nfun),
                    thetag = array(0, dim = c(nbasis + 1, nfun, smcmc)), betag = matrix(0, smcmc, nparw), alphag = matrix(0, smcmc, nfun), sigma2g = numeric(smcmc),
                    psig = array(0, c(smcmc, maxNext, nfun)), omegag = array(0, c(smcmc, maxNext, nfun)), nug = matrix(0, smcmc, nobs),
                    fxgridg = array(0, dim = c(nint + 1, nfun, smcmc)), fxobsg = array(0, dim = c(nobs, nfun, smcmc)), yestg = matrix(0, smcmc, nobs),
                    wbg = matrix(0, smcmc, nobs), loglikeg = numeric(smcmc), logpriorg = numeric(smcmc), imodmetg = as.integer(numeric(1)),
                    pmetg = numeric(nfun), NAOK = TRUE, PACKAGE = "bsamGP")
  })

  mcmc.draws <- list()
  mcmc.draws$zeta <- foo$zetag
  mcmc.draws$tau2 <- foo$tau2g
  mcmc.draws$tau <- apply(foo$tau2g, 2, sqrt)
  mcmc.draws$alpha <- foo$alphag
  mcmc.draws$psi <- foo$psig
  mcmc.draws$omega <- foo$omegag
  mcmc.draws$gamma <- foo$gammag
  mcmc.draws$lngamma <- apply(foo$gammag, 2, log)
  mcmc.draws$theta <- foo$thetag
  mcmc.draws$beta <- foo$betag
  mcmc.draws$sigma2 <- foo$sigma2g
  mcmc.draws$sigma <- sqrt(foo$sigma2g)
  mcmc.draws$nu <- foo$nug

  loglik.draws <- list()
  loglik.draws$loglike <- foo$loglikeg
  loglik.draws$logprior <- foo$logpriorg
  loglik.draws$logjoint <- foo$loglikeg + foo$logpriorg

  fit.draws <- list()
  fit.draws$xgrid <- xgrid
  fit.draws$fxobs <- foo$fxobsg
  fit.draws$fxgrid <- foo$fxgridg
  fit.draws$yhat <- foo$yestg
  fit.draws$wbeta <- foo$wbg

  post.est <- list()
  betam <- apply(foo$betag, 2, mean)
  betas <- apply(foo$betag, 2, sd)
  post.est$betam <- betam
  post.est$betas <- betas

  sigma2m <- mean(foo$sigma2g)
  sigma2s <- sd(foo$sigma2g)
  post.est$sigma2m <- sigma2m
  post.est$sigma2s <- sigma2s

  sigmag <- sqrt(foo$sigma2g)
  sigmam <- mean(sigmag)
  sigmas <- sd(sigmag)
  post.est$sigmam <- sigmam
  post.est$sigmas <- sigmas

  thetam <- apply(foo$thetag, c(1, 2), mean)
  thetas <- apply(foo$thetag, c(1, 2), sd)
  post.est$thetam <- thetam
  post.est$thetas <- thetas

  tau2m <- colMeans(foo$tau2g)
  tau2s <- apply(foo$tau2g, 2, sd)
  post.est$tau2m <- tau2m
  post.est$tau2s <- tau2s

  taug <- sqrt(foo$tau2g)
  taum <- colMeans(taug)
  taus <- apply(taug, 2, sd)
  post.est$taum <- taum
  post.est$taus <- taus

  gammam <- colMeans(foo$gammag)
  gammas <- apply(foo$gammag, 2, sd)
  post.est$gammam <- gammam
  post.est$gammas <- gammas

  alpham <- colMeans(foo$alphag)
  alphas <- apply(foo$alphag, 2, sd)
  post.est$alpham <- alpham
  post.est$alphas <- alphas

  psim <- apply(foo$psig, c(2,3), mean)
  psis <- apply(foo$psig, c(2,3), sd)
  post.est$psim <- psim
  post.est$psis <- psis

  omegam <- apply(foo$omegag, c(2,3), mean)
  omegas <- apply(foo$omegag, c(2,3), sd)
  post.est$omegam <- omegam
  post.est$omegas <- omegas

  zetam <- colMeans(foo$zetag)
  zetas <- apply(foo$zetag, 2, sd)
  post.est$zetam <- zetam
  post.est$zetas <- zetas

  if (marginal.likelihood) {
    sfact <- 0.75
    vfact <- sfact^2

    beta_mn <- betam
    beta_sn <- sfact * betas
    beta_cov <- matrix(beta_sn, nparw, nparw) * (cor(foo$betag)) * matrix(beta_sn, nparw, nparw, byrow = TRUE)
    beta_covi <- solve(beta_cov)
    lndetbcov <- log(det(beta_cov))

    sigma2_rn <- 2 * (2 + (sigma2m/(sfact * sigma2s))^2)
    sigma2_sn <- sigma2m * (sigma2_rn - 2)

    theta_mn <- thetam
    theta_sn <- sfact * thetas
    theta_vn <- theta_sn^2
    theta0_mn <- thetam[1, ]
    theta0_sn <- sfact * thetas[1, ]
    theta0_vn <- theta0_sn^2
    theta0_lnpn <- numeric(nfun)
    for (i in 1:nfun) {
      if (fmodel[i] > 1) {
        theta0_lnpn[i] <- pnorm(-theta0_mn[i]/theta0_sn[i], lower.tail = FALSE, log.p = TRUE)
      }
    }

    tau2_rn <- 2 * (2 + (tau2m/(sfact * tau2s))^2)
    tau2_sn <- tau2m * (tau2_rn - 2)

    gamma_mn <- gammam
    gamma_sn <- (sfact * gammas)
    gamma_vn <- gamma_sn^2
    gamma_lnpn <- pnorm(-gamma_mn/gamma_sn, lower.tail = FALSE, log.p = TRUE)

    alpha_mn <- alpham
    alpha_sn <- sfact * alphas
    alpha_vn <- alpha_sn^2
    alpha_lnpn <- numeric(nfun)
    for (i in 1:nfun) {
      if (fmodel[i] > 2) {
        if (fpm[i] == 1) {
          alpha_lnpn[i] <- pnorm(-alpha_mn[i]/alpha_sn[i], lower.tail = FALSE, log.p = TRUE)
        } else {
          alpha_lnpn[i] <- pnorm(-alpha_mn[i]/alpha_sn[i], lower.tail = FALSE, log.p = TRUE)
        }
      }
    }

    psi_mn <- psim
    psi_sn <- psis * sfact
    psi_vn <- psi_sn^2

    omega_mn <- omegam
    omega_sn <- sfact * omegas
    omega_vn <- omega_sn^2

    psi_lnpn <- numeric(nfun)
    omega_lnpn <- numeric(nfun)
    for (i in 1:nfun) {
      if ((4 < fmodel[i]) && (fmodel[i] < 8)) {
        if (iflagpsi) {
          if (psi_sn[1,i] <= 0) {
            psi_sn[1,i] <- 0.001
            psi_vn[1,i] <- psi_sn[1,i]^2
          }
          psi_lnpn[i] <- log(1 - pnorm(-psi_mn[1,i]/psi_sn[1,i]))
        }
        if (omega_sn[1,i] <= 0) {
          omega_sn[1,i] <- 0.01
          omega_vn[1,i] <- omega_sn[1,i]^2
        }
        omega_lnpn[i] <- log(pnorm((xmax[i] - omegam[1,i])/(omega_sn[1,i])) - pnorm((xmin[i] - omegam[1,i])/(omega_sn[1,i])))
      }
      if (fmodel[i] == 8) {
        if(iflagpsi) {
          psi_sn[psi_sn[,i] <= 0, i] <- 0.001
          psi_vn[,i] <- psi_sn[,i]^2
          psi_lnpn[i] <- sum(log(1 - pnorm(-psi_mn[,i]/psi_sn[,i])))
        }
        omega_sn[omega_sn[,i] <= 0, i] <- 0.01
        omega_vn[,i] <- omega_sn[,i]^2

        if(nExtremes[i]<3) {
          omega_lnpn[i] <- sum( log(pnorm( (omegam[2,i]-omegam[1,i])/(omega_sn[1,i])) - pnorm( (xmin[i]-omegam[1,i])/(omega_sn[1,i]))),
                  log(pnorm( (xmax[i]-omegam[nExtremes[i],i])/(omega_sn[nExtremes[i],i])) - pnorm( (omegam[nExtremes[i]-1,i]-omegam[nExtremes[i],i])/(omega_sn[nExtremes[i],i]))) )
        } else {
          omega_lnpn[i] <- sum( log(pnorm( (omegam[2,i]-omegam[1,i])/(omega_sn[1,i])) - pnorm( (xmin[i]-omegam[1,i])/(omega_sn[1,i]))),
              sapply(2:(nExtremes[i]-1), function(k) log(pnorm( (omegam[k+1,i]-omegam[k,i])/(omega_sn[k,i])) - pnorm( (omegam[k-1,i]-omegam[k,i])/(omega_sn[k,i])))),
              log(pnorm( (xmax[i]-omegam[nExtremes[i],i])/(omega_sn[nExtremes[i],i])) - pnorm( (omegam[nExtremes[i]-1,i]-omegam[nExtremes[i],i])/(omega_sn[nExtremes[i],i]))) )
        }
      }
    }

    logg <- .Fortran("bsaramgetlogg", as.integer(fmodel), as.matrix(mcmc.draws$beta), as.double(mcmc.draws$sigma2), as.array(mcmc.draws$theta),
                     as.matrix(mcmc.draws$tau2), as.matrix(mcmc.draws$gamma), as.matrix(mcmc.draws$alpha), as.array(mcmc.draws$psi),
                     as.array(mcmc.draws$omega), as.integer(smcmc), as.integer(nparw), as.integer(nfun), as.integer(nbasis), as.integer(nExtremes), as.integer(maxNext), as.integer(iflagpsi),
                     as.double(beta_mn), as.matrix(beta_covi), as.double(lndetbcov), as.double(sigma2_rn), as.double(sigma2_sn), as.double(theta_mn),
                     as.double(theta_sn), as.double(theta0_lnpn), as.double(tau2_rn), as.double(tau2_sn), as.double(gamma_mn), as.double(gamma_vn),
                     as.double(gamma_lnpn), as.double(alpha_mn), as.double(alpha_vn), as.double(alpha_lnpn), as.matrix(psi_mn), as.matrix(psi_vn),
                     as.double(psi_lnpn), as.matrix(omega_mn), as.matrix(omega_vn), as.double(omega_lnpn), logg = numeric(smcmc),
                     NAOK = TRUE, PACKAGE = "bsamGP")$logg
    loglik.draws$logg <- logg

    logjointg <- foo$loglikeg + foo$logpriorg
    ratiog <- logg - logjointg
    mratio <- max(ratiog)
    lilg <- exp(ratiog - mratio)
    lilm <- -mratio - log(mean(lilg))

    a <- max(-foo$loglikeg)
    lilnr <- -a - log(mean(exp(-foo$loglikeg - a)))

    if (spm.adequacy) {
      fout.lm <- .Fortran("bqreg", as.double(yobs), as.matrix(wdata), as.double(p), as.double(betam), as.double(beta_m0),
                          as.matrix(beta_v0), as.double(sigma2m), as.double(sigma2_m0), as.double(sigma2_v0), as.integer(nobs), as.integer(nparw),
                          as.integer(nblow), as.integer(nskip), as.integer(smcmc), betaps = matrix(0, smcmc, nparw), sigmasqps = numeric(smcmc),
                          loglikeps = numeric(smcmc), logpriorps = numeric(smcmc), NAOK = TRUE, PACKAGE = "bsamGP")
      logjoint.lm <- fout.lm$logpriorps + fout.lm$loglikeps

      beta_mn <- colMeans(fout.lm$betaps)
      beta_sn <- sfact * apply(fout.lm$betaps, 2, sd)
      beta_cov <- matrix(beta_sn, nparw, nparw) * (cor(fout.lm$betaps)) * matrix(beta_sn, nparw, nparw)
      beta_covi <- solve(beta_cov)
      lndetbcov <- log(det(beta_cov))

      sigma2_rn <- 2 * (2 + (mean(fout.lm$sigmasqps)/(sfact * sd(fout.lm$sigmasqps)))^2)
      sigma2_sn <- mean(fout.lm$sigmasqps) * (sigma2_rn - 2)

      logg.lm <- .Fortran("bqreggetlogg", as.matrix(fout.lm$betaps), as.double(fout.lm$sigmasqps), as.integer(smcmc),
                          as.integer(nparw), as.double(beta_mn), as.matrix(beta_covi), as.double(lndetbcov), as.double(sigma2_rn), as.double(sigma2_sn),
                          logg = numeric(smcmc), NAOK = TRUE, PACKAGE = "bsamGP")$logg

      ratiog.lm <- logg.lm - logjoint.lm
      mratio.lm <- max(ratiog.lm)
      lilg.lm <- exp(ratiog.lm - mratio.lm)
      lil.lm <- -mratio.lm - log(mean(lilg.lm))
    }
  }

  ym <- colMeans(foo$yestg)
  rsquarey <- cor(cbind(yobs, ym))^2
  rsquarey <- rsquarey[1, 2]

  res.out <- list()
  res.out$call <- cl
  res.out$model <- "bsaq"
  res.out$y <- yobs
  res.out$w <- wdata
  res.out$x <- xobs
  res.out$p <- p
  res.out$xmin <- xmin
  res.out$xmax <- xmax
  res.out$n <- nobs
  res.out$ndimw <- ndimw
  res.out$nparw <- nparw
  res.out$nint <- nint
  res.out$nbasis <- nbasis

  res.out$yname <- yname
  res.out$wnames <- wnames
  res.out$xname <- xname

  res.out$shape <- shape
  res.out$fshape <- fshape
  res.out$fmodel <- fmodel
  res.out$fpm <- fpm
  res.out$nfun <- nfun
  res.out$nExtremes <- nExtremes
  
  res.out$prior <- privals

  res.out$mcmctime <- mcmctime
  res.out$mcmc <- mcvals
  res.out$pmet <- foo$pmetg
  res.out$imodmet <- foo$imodmetg

  res.out$mcmc.draws <- mcmc.draws
  res.out$fit.draws <- fit.draws
  res.out$loglik.draws <- loglik.draws

  res.out$post.est <- post.est

  res.out$marglik <- marginal.likelihood
  res.out$spmadeq <- spm.adequacy
  if (marginal.likelihood) {
    if (spm.adequacy)
      res.out$lmarg.lm <- lil.lm[1]
    res.out$lmarg.gd <- lilm
    res.out$lmarg.nr <- lilnr
  }

  res.out$rsquarey <- rsquarey

  class(res.out) <- "bsam"
  res.out
}
