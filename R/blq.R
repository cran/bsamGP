"blq" <- function(y, w, p, mcmc = list(), prior = list(), marginal.likelihood = TRUE) {
  cl <- match.call()
  yobs <- y
  if (missing(w)) {
    wdata <- NULL
  } else {
    wdata <- w
  }

  if (missing(p)) {
    p <- 0.5
  }

  if (!is.matrix(yobs))
    yobs <- as.matrix(yobs)
  yname <- colnames(yobs)
  if (is.null(yname))
    yname <- "y"
  colnames(yobs) <- yname
  nobs <- nrow(yobs)

  if (!is.matrix(wdata))
    wdata <- as.matrix(wdata)
  wnames <- colnames(wdata)
  if (is.null(wnames))
    wnames <- paste("w", 1:ncol(wdata), sep = "")
  wdata <- cbind(1, wdata)
  wnames <- c("const", wnames)
  colnames(wdata) <- wnames
  ndimw <- ncol(wdata) - 1
  nparw <- ndimw + 1

  privals <- list(beta_m0 = numeric(nparw), beta_v0 = diag(100, nparw),
                  sigma2_m0 = 1, sigma2_v0 = 1000)
  privals[names(prior)] <- prior

  beta_m0 <- privals$beta_m0
  beta_v0 <- privals$beta_v0

  sigma2_m0 <- privals$sigma2_m0
  sigma2_v0 <- privals$sigma2_v0

  mcvals <- list(nblow = 1000, smcmc = 1000, nskip = 1)
  mcvals[names(mcmc)] <- mcmc
  nblow <- mcvals$nblow
  smcmc <- mcvals$smcmc
  nskip <- mcvals$nskip

  mcmctime <- system.time({
    fout <- .Fortran("bqreg", as.double(yobs), as.matrix(wdata), as.double(p),
                     as.double(beta_m0), as.double(beta_m0), as.matrix(beta_v0),
                     as.double(sigma2_m0), as.double(sigma2_m0), as.double(sigma2_v0),
                     as.integer(nobs), as.integer(nparw),
                     as.integer(nblow), as.integer(nskip), as.integer(smcmc),
                     betaps = matrix(0, smcmc, nparw), sigmasqps = numeric(smcmc),
                     loglikeps = numeric(smcmc), logpriorps = numeric(smcmc),
                     NAOK = TRUE, PACKAGE = "bsamGP")
  })

  mcmc.draws <- list()
  mcmc.draws$beta <- fout$betaps
  mcmc.draws$sigma2 <- fout$sigmasqps
  mcmc.draws$sigma <- sqrt(fout$sigmasqps)

  loglik.draws <- list()
  loglik.draws$loglike <- fout$loglikeps
  loglik.draws$logprior <- fout$logpriorps
  loglik.draws$logjoint <- fout$loglikeps + fout$logpriorps

  fit.draws <- list()
  fit.draws$wbeta <- tcrossprod(fout$betaps, wdata)

  post.est <- list()
  betam <- apply(fout$betaps, 2, mean)
  betas <- apply(fout$betaps, 2, sd)
  post.est$betam <- betam
  post.est$betas <- betas

  sigma2m <- mean(fout$sigmasqps)
  sigma2s <- sd(fout$sigmasqps)
  post.est$sigma2m <- sigma2m
  post.est$sigma2s <- sigma2s

  sigmag <- sqrt(fout$sigmasqps)
  sigmam <- mean(sigmag)
  sigmas <- sd(sigmag)
  post.est$sigmam <- sigmam
  post.est$sigmas <- sigmas

  if (marginal.likelihood) {
    sfact <- 0.75
    vfact <- sfact^2

    beta_mn <- colMeans(fout$betaps)
    beta_sn <- sfact * apply(fout$betaps, 2, sd)
    beta_cov <- matrix(beta_sn, nparw, nparw) * (cor(fout$betaps)) * matrix(beta_sn, nparw, nparw)
    beta_covi <- solve(beta_cov)
    lndetbcov <- log(det(beta_cov))

    sigma2_rn <- 2 * (2 + (mean(fout$sigmasqps)/(sfact * sd(fout$sigmasqps)))^2)
    sigma2_sn <- mean(fout$sigmasqps) * (sigma2_rn - 2)

    logg.lm <- .Fortran("bqreggetlogg", as.matrix(fout$betaps), as.double(fout$sigmasqps),
                        as.integer(smcmc), as.integer(nparw),
                        as.double(beta_mn), as.matrix(beta_covi), as.double(lndetbcov),
                        as.double(sigma2_rn), as.double(sigma2_sn), logg = numeric(smcmc),
                        NAOK = TRUE, PACKAGE = "bsamGP")$logg
    logjoint.lm <- fout$logpriorps + fout$loglikeps
    ratiog.lm <- logg.lm - logjoint.lm
    mratio.lm <- max(ratiog.lm)
    lilg.lm <- exp(ratiog.lm - mratio.lm)
    lil <- -mratio.lm - log(mean(lilg.lm))
  }

  ym <- colMeans(tcrossprod(fout$betaps, wdata))
  rsquarey <- cor(cbind(yobs, ym))^2
  rsquarey <- rsquarey[1, 2]

  res.out <- list()
  res.out$call <- cl
  res.out$model <- "blq"
  res.out$y <- yobs
  res.out$w <- wdata
  res.out$p <- p
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
  res.out$loglik.draws <- loglik.draws

  res.out$post.est <- post.est

  res.out$marglik <- marginal.likelihood
  if (marginal.likelihood) {
    res.out$lmarg <- lil[1]
  }

  res.out$rsquarey <- rsquarey

  class(res.out) <- "blm"
  res.out
}
