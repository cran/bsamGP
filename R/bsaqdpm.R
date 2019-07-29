"bsaqdpm" <- function(formula, xmin, xmax, p, nbasis, nint, mcmc = list(), prior = list(), egrid, ngrid = 500,
                      shape = c("Free", "Increasing", "Decreasing", "IncreasingConvex", "DecreasingConcave",
                                "IncreasingConcave", "DecreasingConvex", "IncreasingS", "DecreasingS",
                                "IncreasingRotatedS", "DecreasingRotatedS", "InvertedU", "Ushape"),
                      verbose = FALSE) {
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

  if (missing(egrid)) {
    egrid <- seq(-10, 10, length = ngrid)
  } else {
    ngrid <- length(egrid)
  }

  privals <- list(iflagprior = 0, theta0_m0 = 0, theta0_s0 = 100, tau2_m0 = 1, tau2_v0 = 100, w0 = 2, beta_m0 = numeric(nparw),
                  beta_v0 = diag(100, nparw), sigma2_r0 = 4, sigma2_s0 = 0.04 * range(yobs)^2, alpha_m0 = 3, alpha_s0 = 50, iflagpsi = 1,
                  psifixed = 100, omega_m0 = (xmin + xmax)/2, omega_s0 = xrange/8, tmass_a = 2, tmass_b = 4)
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

  sigma2_r0 <- privals$sigma2_r0
  sigma2_s0 <- privals$sigma2_s0
  tmass_a <- privals$tmass_a
  tmass_b <- privals$tmass_b

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
    foo <- .Fortran("bsaqamdpscale", as.integer(verbose), as.double(yobs), as.matrix(wdata), as.matrix(xobs), as.double(egrid), as.integer(nobs),
                    as.integer(ngrid), as.integer(nparw), as.integer(nfun), as.integer(nbasis), as.integer(nint), as.integer(fmodel),
                    as.double(fpm), as.double(p), as.double(theta0_m0), as.double(theta0_s0), as.double(tau2_m0), as.double(tau2_v0),
                    as.double(w0), as.double(beta_m0), as.matrix(beta_v0), as.double(alpha_m0), as.double(alpha_s0), as.double(psi_m0),
                    as.double(psi_s0), as.double(psifixed), as.double(omega_m0), as.double(omega_s0), as.double(sigma2_r0), as.double(sigma2_s0),
                    as.double(tmass_a), as.double(tmass_b), as.integer(iflagprior), as.integer(iflagpsi), as.integer(maxmodmet), as.integer(nblow0),
                    as.integer(nblow), as.integer(smcmc), as.integer(nskip), as.integer(ndisp), zetag = matrix(0, smcmc, nfun),
                    tau2g = matrix(0, smcmc, nfun), gammag = matrix(0, smcmc, nfun), thetag = array(0, dim = c(nbasis + 1, nfun, smcmc)),
                    betag = matrix(0, smcmc, nparw), alphag = matrix(0, smcmc, nfun), sigma2g = matrix(0, nrow = smcmc, ncol = nobs),
                    psig = matrix(0, smcmc, nfun), omegag = matrix(0, smcmc, nfun), nug = matrix(0, smcmc, nobs),
                    fxgridg = array(0, dim = c(nint + 1, nfun, smcmc)), fxobsg = array(0, dim = c(nobs, nfun, smcmc)),
                    yestg = matrix(0, smcmc, nobs), wbg = matrix(0, smcmc, nobs), tmassg = numeric(smcmc), config = as.integer(matrix(0, smcmc, nobs)),
                    nclassg = as.integer(numeric(smcmc)), edensg = matrix(0, smcmc, ngrid), invlikeg = matrix(0, nrow = smcmc, ncol = nobs),
                    imodmetg = as.integer(numeric(1)), pmetg = numeric(nfun), NAOK = TRUE, PACKAGE = "bsamGP")
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
  mcmc.draws$nu <- foo$nug

  dpm.draws <- list()
  dpm.draws$tmass <- foo$tmassg
  dpm.draws$config <- foo$config
  dpm.draws$nclass <- foo$nclassg
  dpm.draws$edens <- foo$edensg
  dpm.draws$egrid <- egrid
  dpm.draws$ngrid <- ngrid

  lik.draws <- list()
  lik.draws$invlike <- foo$invlikeg

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

  sigma2m <- colMeans(foo$sigma2g)
  sigma2s <- apply(foo$sigma2g, 2, sd)
  post.est$sigma2m <- sigma2m
  post.est$sigma2s <- sigma2s

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

  ym <- colMeans(foo$yestg)
  rsquarey <- cor(cbind(yobs, ym))^2
  rsquarey <- rsquarey[1, 2]

  res.out <- list()
  res.out$call <- cl
  res.out$model <- "bsaqdpm"
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
  res.out$location <- FALSE

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
  res.out$lik.draws <- lik.draws
  res.out$dpm.draws <- dpm.draws

  res.out$post.est <- post.est

  res.out$rsquarey <- rsquarey
  res.out$cpo <- 1/colMeans(foo$invlikeg)
  res.out$lpml <- sum(log(res.out$cpo))

  class(res.out) <- "bsamdpm"
  res.out
}
