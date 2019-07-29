"blr" <- function(formula, data = NULL, mcmc = list(), prior = list(), marginal.likelihood = TRUE) {
  cl <- match.call()
  
  ywdata <- parse.formula(formula, data = data)
  yobs <- ywdata[[1]]
  yname <- ywdata[[2]]
  wdata <- ywdata[[3]]
  wnames <- ywdata[[4]]

  nobs <- nrow(yobs)
  nparw <- ncol(wdata)
  ndimw <- nparw - 1

  privals <- list(beta_m0 = numeric(nparw), beta_v0 = diag(100, nparw),
                  sigma2_m0 = 1, sigma2_v0 = 1000)
  privals[names(prior)] <- prior

  beta_m0 <- privals$beta_m0
  beta_v0 <- privals$beta_v0
  beta_iv0 <- solve(beta_v0)
  beta_iv0m0 <- beta_iv0 %*% beta_m0

  sigma2_m0 <- privals$sigma2_m0
  sigma2_v0 <- privals$sigma2_v0
  sigma2_r0 <- 2 * (2 + sigma2_m0^2/sigma2_v0)
  sigma2_s0 <- sigma2_m0 * (sigma2_r0 - 2)

  mcvals <- list(nblow = 1000, smcmc = 1000, nskip = 1)
  mcvals[names(mcmc)] <- mcmc
  nblow <- mcvals$nblow
  smcmc <- mcvals$smcmc
  nskip <- mcvals$nskip
  nmcmc <- nblow + nskip * smcmc
  s.ind <- seq(nblow + 1, nmcmc, nskip)

  mcmctime <- system.time({
    wtw <- crossprod(wdata)
    wty <- crossprod(wdata, yobs)
    yty <- crossprod(yobs)

    beta_ivn <- wtw + beta_iv0
    beta_vn <- solve(beta_ivn)
    beta_mn <- beta_vn %*% (wty + beta_iv0m0)
    beta_ivnmn <- beta_ivn %*% beta_mn

    sigma2_rn <- sigma2_r0 + nobs
    sigma2_sn <- sigma2_s0 + yty + crossprod(beta_m0, beta_iv0m0) -
      crossprod(beta_mn, beta_ivnmn)

    sigma2g <- 1/rgamma(nmcmc, sigma2_rn/2, sigma2_sn/2)

    betag <- .Fortran("blreg", as.double(sigma2g), as.double(beta_mn), as.matrix(beta_vn),
                      as.integer(nparw), as.integer(nmcmc),
                      betag = matrix(0, nmcmc, nparw), NAOK = TRUE, PACKAGE = "bsamGP")$betag

    sigma2g <- sigma2g[s.ind]
    betag <- betag[s.ind, ]
  })

  mcmc.draws <- list()
  mcmc.draws$beta <- betag
  mcmc.draws$sigma2 <- sigma2g
  mcmc.draws$sigma <- sqrt(sigma2g)

  fit.draws <- list()
  fit.draws$wbeta <- tcrossprod(betag, wdata)

  post.est <- list()
  betam <- apply(betag, 2, mean)
  betas <- apply(betag, 2, sd)
  post.est$betam <- betam
  post.est$betas <- betas

  sigma2m <- mean(sigma2g)
  sigma2s <- sd(sigma2g)
  post.est$sigma2m <- sigma2m
  post.est$sigma2s <- sigma2s

  sigmag <- sqrt(sigma2g)
  sigmam <- mean(sigmag)
  sigmas <- sd(sigmag)
  post.est$sigmam <- sigmam
  post.est$sigmas <- sigmas

  if (marginal.likelihood) {
    a0 <- sigma2_r0/2
    b0 <- sigma2_s0/2
    an <- sigma2_rn/2
    bn <- sigma2_sn/2
    lil <- -nobs * log(2 * pi)/2 + log(det(beta_iv0))/2 - log(det(beta_ivn))/2 +
      a0 * log(b0) - an * log(bn) + lgamma(an) - lgamma(a0)
  }

  ym <- colMeans(fit.draws$wbeta)
  rsquarey <- cor(cbind(yobs, ym))^2
  rsquarey <- rsquarey[1, 2]

  res.out <- list()
  res.out$call <- cl
  res.out$model <- "blr"
  res.out$y <- yobs
  res.out$w <- wdata
  res.out$n <- nobs
  res.out$ndimw <- ndimw
  res.out$nparw <- nparw
  res.out$yname <- yname
  res.out$wnames <- wnames
  res.out$prior <- privals
  res.out$mcmctime <- mcmctime
  res.out$mcmc <- mcvals
  res.out$mcmc.draws <- mcmc.draws
  res.out$fit.draws <- fit.draws
  res.out$post.est <- post.est
  res.out$marglik <- marginal.likelihood
  if (marginal.likelihood) {
    res.out$lmarg <- lil[1]
  }
  res.out$rsquarey <- rsquarey

  class(res.out) <- "blm"
  res.out
}
