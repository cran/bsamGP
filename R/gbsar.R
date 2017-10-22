"gbsar" <- function(formula, xmin, xmax, family, link, nbasis, nint, mcmc = list(), prior = list(),
                    shape = c("Free", "Increasing", "Decreasing", "IncreasingConvex", "DecreasingConcave",
                              "IncreasingConcave", "DecreasingConvex", "IncreasingS", "DecreasingS",
                              "IncreasingRotatedS", "DecreasingRotatedS", "InvertedU", "Ushape"),
                    marginal.likelihood = TRUE, algorithm = c("AM", "KS")) {
  cl <- match.call()

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

  if (family %in% c("bernoulli", "poisson", "negative.binomial", "poisson.gamma")) {
    if (link %in% c("probit", "logit", "log")) {
      if (family %in% c("poisson", "negative.binomial", "poisson.gamma") && link != "log") {
        stop(paste(link, " link is not recognized in ", family, " model", sep = ""))
      }
      if (family == "bernoulli" && link == "log") {
        stop(paste(link, " link is not recognized in ", family, " model", sep = ""))
      }
    } else {
      stop(paste(link, " link is not recognized.", sep = ""))
    }
  } else {
    stop(paste(family, " model is not recognized.", sep = ""))
  }

  if (family == "bernoulli") {
    chk <- unique(yobs)
    if (length(chk) != 2L) {
      stop("bernoulli models only allow binary responses")
    } else {
      yobs <- factor(yobs, labels = 0:1)
      yobs <- as.integer(paste(yobs))
    }
  }
  if (family %in% c("poisson", "negative.binomial", "poisson.gamma")) {
    chk <- (abs(yobs - round(yobs)) < .Machine$double.eps^0.5)
    if (sum(chk) != length(yobs))
      stop(paste(family, " model only allow count responses.", sep = ""))
  }
  if (family == "bernoulli" && link == "logit") {
    algorithm = match.arg(algorithm)
  }

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

  privals <- list(iflagprior = 0, theta0_m0 = 0, theta0_s0 = 100, tau2_m0 = 1, tau2_v0 = 100, w0 = 2, beta_m0 = numeric(nparw),
                  beta_v0 = diag(100, nparw), alpha_m0 = 3, alpha_s0 = 50, iflagpsi = 1, psifixed = 100, omega_m0 = (xmin + xmax)/2,
                  omega_s0 = xrange/8, kappa_m0 = 1, kappa_v0 = 100)
  privals[names(prior)] <- prior

  iflagprior <- privals$iflagprior
  theta0_m0 <- privals$theta0_m0
  theta0_s0 <- privals$theta0_s0
  tau2_m0 <- privals$tau2_m0
  tau2_v0 <- privals$tau2_v0
  w0 <- privals$w0
  beta_m0 <- privals$beta_m0
  beta_v0 <- privals$beta_v0
  alpha_m0 <- privals$alpha_m0
  alpha_s0 <- privals$alpha_s0
  iflagpsi <- privals$iflagpsi
  psifixed <- privals$psifixed
  psi_m0 <- psifixed
  psi_s0 <- 10 * psifixed
  omega_m0 <- privals$omega_m0
  omega_s0 <- privals$omega_s0
  kappa_m0 <- privals$kappa_m0
  kappa_v0 <- privals$kappa_v0

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

  stime <- proc.time()
  if (family == "bernoulli" && link == "probit") {
    foo <- .Fortran("gbsarprobitAC", as.integer(yobs), as.matrix(wdata), as.matrix(xobs), as.integer(nobs), as.integer(nparw),
                    as.integer(nfun), as.integer(nbasis), as.integer(nint), as.integer(fmodel), as.double(fpm), as.double(theta0_m0),
                    as.double(theta0_s0), as.double(tau2_m0), as.double(tau2_v0), as.double(w0), as.double(beta_m0), as.matrix(beta_v0),
                    as.double(alpha_m0), as.double(alpha_s0), as.double(psi_m0), as.double(psi_s0), as.double(psifixed), as.double(omega_m0),
                    as.double(omega_s0), as.integer(iflagprior), as.integer(iflagpsi), as.integer(maxmodmet), as.integer(nblow0), as.integer(nblow),
                    as.integer(smcmc), as.integer(nskip), as.integer(ndisp), zetag = matrix(0, smcmc, nfun), tau2g = matrix(0, smcmc, nfun),
                    gammag = matrix(0, smcmc, nfun), thetag = array(0, dim = c(nbasis + 1, nfun, smcmc)), betag = matrix(0, smcmc, nparw),
                    alphag = matrix(0, smcmc, nfun), psig = matrix(0, smcmc, nfun), omegag = matrix(0, smcmc, nfun),
                    fxgridg = array(0, dim = c(nint + 1, nfun, smcmc)), fxobsg = array(0, dim = c(nobs, nfun, smcmc)), muhatg = matrix(0, smcmc, nobs),
                    wbg = matrix(0, smcmc, nobs), loglikeg = numeric(smcmc), logpriorg = numeric(smcmc), imodmetg = as.integer(numeric(1)),
                    pmetg = numeric(nfun), NAOK = TRUE, PACKAGE = "bsamGP")
  }
  if (family == "bernoulli" && link == "logit") {
    if (algorithm == "KS") {
      foo <- .Fortran("gbsarlogitKS", as.integer(yobs), as.matrix(wdata), as.matrix(xobs), as.integer(nobs), as.integer(nparw),
                      as.integer(nfun), as.integer(nbasis), as.integer(nint), as.integer(fmodel), as.double(fpm), as.double(theta0_m0),
                      as.double(theta0_s0), as.double(tau2_m0), as.double(tau2_v0), as.double(w0), as.double(beta_m0), as.matrix(beta_v0),
                      as.double(alpha_m0), as.double(alpha_s0), as.double(psi_m0), as.double(psi_s0), as.double(psifixed), as.double(omega_m0),
                      as.double(omega_s0), as.integer(iflagprior), as.integer(iflagpsi), as.integer(maxmodmet), as.integer(nblow0), as.integer(nblow),
                      as.integer(smcmc), as.integer(nskip), as.integer(ndisp), zetag = matrix(0, smcmc, nfun), tau2g = matrix(0, smcmc, nfun),
                      gammag = matrix(0, smcmc, nfun), thetag = array(0, dim = c(nbasis + 1, nfun, smcmc)), betag = matrix(0, smcmc, nparw),
                      alphag = matrix(0, smcmc, nfun), psig = matrix(0, smcmc, nfun), omegag = matrix(0, smcmc, nfun),
                      fxgridg = array(0, dim = c(nint + 1, nfun, smcmc)), fxobsg = array(0, dim = c(nobs, nfun, smcmc)), muhatg = matrix(0, smcmc, nobs),
                      wbg = matrix(0, smcmc, nobs), loglikeg = numeric(smcmc), logpriorg = numeric(smcmc), imodmetg = as.integer(numeric(1)),
                      pmetg = numeric(nfun), NAOK = TRUE, PACKAGE = "bsamGP")
    } else {
      foo <- .Fortran("gbsarlogitMH", as.integer(yobs), as.matrix(wdata), as.matrix(xobs), as.integer(nobs), as.integer(nparw),
                      as.integer(nfun), as.integer(nbasis), as.integer(nint), as.integer(fmodel), as.double(fpm), as.double(theta0_m0),
                      as.double(theta0_s0), as.double(tau2_m0), as.double(tau2_v0), as.double(w0), as.double(beta_m0), as.matrix(beta_v0),
                      as.double(alpha_m0), as.double(alpha_s0), as.double(psi_m0), as.double(psi_s0), as.double(psifixed), as.double(omega_m0),
                      as.double(omega_s0), as.integer(iflagprior), as.integer(iflagpsi), as.integer(maxmodmet), as.integer(nblow0), as.integer(nblow),
                      as.integer(smcmc), as.integer(nskip), as.integer(ndisp), zetag = matrix(0, smcmc, nfun), tau2g = matrix(0, smcmc, nfun),
                      gammag = matrix(0, smcmc, nfun), thetag = array(0, dim = c(nbasis + 1, nfun, smcmc)), betag = matrix(0, smcmc, nparw),
                      alphag = matrix(0, smcmc, nfun), psig = matrix(0, smcmc, nfun), omegag = matrix(0, smcmc, nfun),
                      fxgridg = array(0, dim = c(nint + 1, nfun, smcmc)), fxobsg = array(0, dim = c(nobs, nfun, smcmc)), muhatg = matrix(0, smcmc, nobs),
                      wbg = matrix(0, smcmc, nobs), loglikeg = numeric(smcmc), logpriorg = numeric(smcmc), imodmetg = as.integer(numeric(1)),
                      pmetg = numeric(nfun), NAOK = TRUE, PACKAGE = "bsamGP")
    }
  }
  if (family == "poisson") {
    foo <- .Fortran("gbsarpoisMH", as.integer(yobs), as.matrix(wdata), as.matrix(xobs), as.integer(nobs), as.integer(nparw),
                    as.integer(nfun), as.integer(nbasis), as.integer(nint), as.integer(fmodel), as.double(fpm), as.double(theta0_m0),
                    as.double(theta0_s0), as.double(tau2_m0), as.double(tau2_v0), as.double(w0), as.double(beta_m0), as.matrix(beta_v0),
                    as.double(alpha_m0), as.double(alpha_s0), as.double(psi_m0), as.double(psi_s0), as.double(psifixed), as.double(omega_m0),
                    as.double(omega_s0), as.integer(iflagprior), as.integer(iflagpsi), as.integer(maxmodmet), as.integer(nblow0), as.integer(nblow),
                    as.integer(smcmc), as.integer(nskip), as.integer(ndisp), zetag = matrix(0, smcmc, nfun), tau2g = matrix(0, smcmc, nfun),
                    gammag = matrix(0, smcmc, nfun), thetag = array(0, dim = c(nbasis + 1, nfun, smcmc)), betag = matrix(0, smcmc, nparw),
                    alphag = matrix(0, smcmc, nfun), psig = matrix(0, smcmc, nfun), omegag = matrix(0, smcmc, nfun),
                    fxgridg = array(0, dim = c(nint + 1, nfun, smcmc)), fxobsg = array(0, dim = c(nobs, nfun, smcmc)), muhatg = matrix(0, smcmc, nobs),
                    wbg = matrix(0, smcmc, nobs), loglikeg = numeric(smcmc), logpriorg = numeric(smcmc), imodmetg = as.integer(numeric(1)),
                    pmetg = numeric(nfun), NAOK = TRUE, PACKAGE = "bsamGP")
  }
  if (family == "negative.binomial") {
    foo <- .Fortran("gbsarnegbinMH", as.integer(yobs), as.matrix(wdata), as.matrix(xobs), as.integer(nobs), as.integer(nparw),
                    as.integer(nfun), as.integer(nbasis), as.integer(nint), as.integer(fmodel), as.double(fpm), as.double(theta0_m0),
                    as.double(theta0_s0), as.double(tau2_m0), as.double(tau2_v0), as.double(w0), as.double(beta_m0), as.matrix(beta_v0),
                    as.double(alpha_m0), as.double(alpha_s0), as.double(psi_m0), as.double(psi_s0), as.double(psifixed), as.double(omega_m0),
                    as.double(omega_s0), as.double(kappa_m0), as.double(kappa_v0), as.integer(iflagprior), as.integer(iflagpsi), as.integer(maxmodmet),
                    as.integer(nblow0), as.integer(nblow), as.integer(smcmc), as.integer(nskip), as.integer(ndisp), zetag = matrix(0, smcmc, nfun),
                    tau2g = matrix(0, smcmc, nfun), gammag = matrix(0, smcmc, nfun), thetag = array(0, dim = c(nbasis + 1, nfun, smcmc)),
                    betag = matrix(0, smcmc, nparw), alphag = matrix(0, smcmc, nfun), psig = matrix(0, smcmc, nfun),
                    omegag = matrix(0, smcmc, nfun), kappag = numeric(smcmc), fxgridg = array(0, dim = c(nint + 1, nfun, smcmc)),
                    fxobsg = array(0, dim = c(nobs, nfun, smcmc)), muhatg = matrix(0, smcmc, nobs), wbg = matrix(0, smcmc, nobs),
                    loglikeg = numeric(smcmc), logpriorg = numeric(smcmc), imodmetg = as.integer(numeric(1)), pmetg = numeric(nfun),
                    NAOK = TRUE, PACKAGE = "bsamGP")
  }
  if (family == "poisson.gamma") {
    foo <- .Fortran("gbsarpoisgammMH", as.integer(yobs), as.matrix(wdata), as.matrix(xobs), as.integer(nobs), as.integer(nparw),
                    as.integer(nfun), as.integer(nbasis), as.integer(nint), as.integer(fmodel), as.double(fpm), as.double(theta0_m0),
                    as.double(theta0_s0), as.double(tau2_m0), as.double(tau2_v0), as.double(w0), as.double(beta_m0), as.matrix(beta_v0),
                    as.double(alpha_m0), as.double(alpha_s0), as.double(psi_m0), as.double(psi_s0), as.double(psifixed), as.double(omega_m0),
                    as.double(omega_s0), as.double(kappa_m0), as.double(kappa_v0), as.integer(iflagprior), as.integer(iflagpsi), as.integer(maxmodmet),
                    as.integer(nblow0), as.integer(nblow), as.integer(smcmc), as.integer(nskip), as.integer(ndisp), zetag = matrix(0, smcmc, nfun),
                    tau2g = matrix(0, smcmc, nfun), gammag = matrix(0, smcmc, nfun), thetag = array(0, dim = c(nbasis + 1, nfun, smcmc)),
                    betag = matrix(0, smcmc, nparw), alphag = matrix(0, smcmc, nfun), psig = matrix(0, smcmc, nfun),
                    omegag = matrix(0, smcmc, nfun), kappag = numeric(smcmc), fxgridg = array(0, dim = c(nint + 1, nfun, smcmc)),
                    fxobsg = array(0, dim = c(nobs, nfun, smcmc)), muhatg = matrix(0, smcmc, nobs), wbg = matrix(0, smcmc, nobs),
                    loglikeg = numeric(smcmc), logpriorg = numeric(smcmc), imodmetg = as.integer(numeric(1)), pmetg = numeric(nfun),
                    NAOK = TRUE, PACKAGE = "bsamGP")
  }
  mcmctime <- proc.time() - stime

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
  if (family == "negative.binomial" || family == "poisson.gamma")
    mcmc.draws$kappa <- foo$kappag

  loglik.draws <- list()
  loglik.draws$loglike <- foo$loglikeg
  loglik.draws$logprior <- foo$logpriorg
  loglik.draws$logjoint <- foo$loglikeg + foo$logpriorg

  fit.draws <- list()
  fit.draws$xgrid <- xgrid
  fit.draws$fxobs <- foo$fxobsg
  fit.draws$fxgrid <- foo$fxgridg
  fit.draws$wbeta <- foo$wbg
  fit.draws$muhat <- foo$muhatg

  post.est <- list()
  betam <- apply(foo$betag, 2, mean)
  betas <- apply(foo$betag, 2, sd)
  post.est$betam <- betam
  post.est$betas <- betas

  if (family == "negative.binomial" || family == "poisson.gamma") {
    kappam <- mean(foo$kappag)
    kappas <- sd(foo$kappag)
    post.est$kappam <- kappam
    post.est$kappas <- kappas
  }

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

  psim <- colMeans(foo$psig)
  psis <- apply(foo$psig, 2, sd)
  post.est$psim <- psim
  post.est$psis <- psis

  omegam <- colMeans(foo$omegag)
  omegas <- apply(foo$omegag, 2, sd)
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
      if (fmodel[i] > 4.5) {
        if (iflagpsi == 1) {
          if (psi_sn[i] <= 0) {
            psi_sn[i] <- 0.001
            psi_vn[i] <- psi_sn[i]^2
          }
          psi_lnpn[i] <- log(1 - pnorm(-psi_mn[i]/psi_sn[i]))
        }
        if (omega_sn[i] <= 0) {
          omega_sn[i] <- 0.01
          omega_vn[i] <- omega_sn[i]^2
        }
        omega_lnpn[i] <- log(pnorm((xmax[i] - omegam[i])/(omega_sn[i])) - pnorm((xmin[i] - omegam[i])/(omega_sn[i])))
      }
    }

    logg <- .Fortran("gbsaramgetlogg", as.integer(fmodel), as.matrix(mcmc.draws$beta), as.array(mcmc.draws$theta), as.matrix(mcmc.draws$tau2),
                     as.matrix(mcmc.draws$gamma), as.matrix(mcmc.draws$alpha), as.matrix(mcmc.draws$psi), as.matrix(mcmc.draws$omega),
                     as.integer(smcmc), as.integer(nparw), as.integer(nfun), as.integer(nbasis), as.integer(iflagpsi), as.double(beta_mn),
                     as.matrix(beta_covi), as.double(lndetbcov), as.double(theta_mn), as.double(theta_sn), as.double(theta0_lnpn), as.double(tau2_rn),
                     as.double(tau2_sn), as.double(gamma_mn), as.double(gamma_vn), as.double(gamma_lnpn), as.double(alpha_mn), as.double(alpha_vn),
                     as.double(alpha_lnpn), as.double(psi_mn), as.double(psi_vn), as.double(psi_lnpn), as.double(omega_mn), as.double(omega_vn),
                     as.double(omega_lnpn), logg = numeric(smcmc), NAOK = TRUE, PACKAGE = "bsamGP")$logg
    loglik.draws$logg <- logg

    logjointg <- foo$loglikeg + foo$logpriorg
    ratiog <- logg - logjointg
    mratio <- max(ratiog)
    lilg <- exp(ratiog - mratio)
    lilm <- -mratio - log(mean(lilg))

    a <- max(-foo$loglikeg)
    lilnr <- -a - log(mean(exp(-foo$loglikeg - a)))
  }

  res.out <- list()
  res.out$call <- cl
  res.out$model <- "gbsar"
  res.out$family <- family
  res.out$link <- link
  res.out$y <- yobs
  res.out$w <- wdata
  res.out$x <- xobs
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
  if (marginal.likelihood) {
    res.out$lmarg.gd <- lilm
    res.out$lmarg.nr <- lilnr
  }

  if (family == "bernoulli" && link == "logit") {
    res.out$algorithm = algorithm
  }

  class(res.out) <- "bsam"
  res.out
}
