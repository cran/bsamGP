"gblr" <- function(formula, data = NULL, family, link, mcmc = list(), prior = list(),
                   marginal.likelihood = TRUE, algorithm = c("AM", "KS"), verbose = FALSE) {
  cl <- match.call()

  ywdata <- parse.formula(formula, data = data)
  yobs <- ywdata[[1]]
  yname <- ywdata[[2]]
  wdata <- ywdata[[3]]
  wnames <- ywdata[[4]]

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

  privals <- list(beta_m0 = numeric(nparw), beta_v0 = diag(100, nparw),
                  kappa_m0 = 1, kappa_v0 = 100)
  privals[names(prior)] <- prior

  beta_m0 <- privals$beta_m0
  beta_v0 <- privals$beta_v0
  kappa_m0 <- privals$kappa_m0
  kappa_v0 <- privals$kappa_v0

  mcvals <- list(nblow = 10000, smcmc = 1000, nskip = 10, ndisp = 1000)
  mcvals[names(mcmc)] <- mcmc
  nblow <- mcvals$nblow
  smcmc <- mcvals$smcmc
  nskip <- mcvals$nskip
  ndisp <- mcvals$ndisp

  stime <- proc.time()
  options(warn=-1)
  if (family == "bernoulli" && link == "probit") {
    betam <- glm(yobs ~ wdata - 1, family = binomial(link = probit))$coef
    fout <- .Fortran("gbprobitAC", as.integer(verbose), as.integer(yobs), as.matrix(wdata),
                     as.double(betam), as.double(beta_m0), as.matrix(beta_v0),
                     as.integer(nobs), as.integer(nparw), as.integer(nblow),
                     as.integer(nskip), as.integer(smcmc), as.integer(ndisp),
                     betaps = matrix(0, smcmc, nparw), loglikeps = numeric(smcmc),
                     logpriorps = numeric(smcmc), NAOK = TRUE, PACKAGE = "bsamGP")
  }
  if (family == "bernoulli" && link == "logit") {
    betam <- glm(yobs ~ wdata - 1, family = binomial(link = logit))$coef
    if (algorithm == 'KS') {
      fout <- .Fortran("gblogitKS", as.integer(verbose), as.integer(yobs), as.matrix(wdata),
                       as.double(betam), as.double(beta_m0), as.matrix(beta_v0),
                       as.integer(nobs), as.integer(nparw),
                       as.integer(nblow), as.integer(nskip), as.integer(smcmc), as.integer(ndisp),
                       betaps = matrix(0, smcmc, nparw), loglikeps = numeric(smcmc),
                       logpriorps = numeric(smcmc), NAOK = TRUE, PACKAGE = "bsamGP")
    } else {
      fout <- .Fortran("gblogitMH", as.integer(verbose), as.integer(yobs), as.matrix(wdata),
                       as.double(betam), as.double(beta_m0), as.matrix(beta_v0),
                       as.integer(nobs), as.integer(nparw),
                       as.integer(nblow), as.integer(nskip), as.integer(smcmc), as.integer(ndisp),
                       betaps = matrix(0, smcmc, nparw), loglikeps = numeric(smcmc),
                       logpriorps = numeric(smcmc), NAOK = TRUE, PACKAGE = "bsamGP")
    }
  }
  if (family == "poisson") {
    glmfit <- glm(y ~ wdata - 1, family = poisson(link = log))
    betam <- glmfit$coef
    mub <- glmfit$coef
    Sb <- vcov(glmfit)

    fout <- .Fortran("gbpoisMH", as.integer(verbose), as.integer(yobs), as.matrix(wdata),
                     as.double(betam), as.double(beta_m0), as.matrix(beta_v0),
                     as.double(mub), as.matrix(Sb), as.integer(nobs),
                     as.integer(nparw), as.integer(nblow), as.integer(nskip), as.integer(smcmc),
                     as.integer(ndisp), betaps = matrix(0, smcmc, nparw),
                     loglikeps = numeric(smcmc), logpriorps = numeric(smcmc), NAOK = TRUE, PACKAGE = "bsamGP")
  }
  if (family == "negative.binomial") {
    betam <- glm.nb(yobs ~ wdata - 1)$coef
    fout <- .Fortran("gbnegbinMH", as.integer(verbose), as.integer(yobs), as.matrix(wdata),
                     as.double(betam), as.double(beta_m0), as.matrix(beta_v0),
                     as.double(kappa_m0), as.double(kappa_v0),
                     as.integer(nobs), as.integer(nparw), as.double(nblow), as.integer(nskip),
                     as.integer(smcmc), as.integer(ndisp), betaps = matrix(0, smcmc, nparw),
                     kappas = numeric(smcmc), loglikeps = numeric(smcmc),
                     logpriorps = numeric(smcmc), NAOK = TRUE, PACKAGE = "bsamGP")
  }
  if (family == "poisson.gamma") {
    betam <- glm.nb(yobs ~ wdata - 1)$coef
    fout <- .Fortran("gbpoisgammMH", as.integer(verbose), as.integer(yobs), as.matrix(wdata),
                     as.double(betam), as.double(beta_m0), as.matrix(beta_v0),
                     as.double(kappa_m0), as.double(kappa_v0), as.integer(nobs),
                     as.integer(nparw), as.double(nblow), as.integer(nskip),
                     as.integer(smcmc), as.integer(ndisp), betaps = matrix(0, smcmc, nparw),
                     kappas = numeric(smcmc), loglikeps = numeric(smcmc),
                     logpriorps = numeric(smcmc), NAOK = TRUE, PACKAGE = "bsamGP")
  }
  options(warn=0)
  mcmctime <- proc.time() - stime

  mcmc.draws <- list()
  mcmc.draws$beta <- fout$betaps
  if (family == "negative.binomial" || family == "poisson.gamma") {
    mcmc.draws$kappa <- fout$kappas
  }

  loglik.draws <- list()
  loglik.draws$loglike <- fout$loglikeps
  loglik.draws$logprior <- fout$logpriorps
  loglik.draws$logjoint <- fout$loglikeps + fout$logpriorps

  fit.draws <- list()
  fit.draws$wbeta <- fout$betaps %*% t(wdata)

  post.est <- list()
  betam <- apply(fout$betaps, 2, mean)
  betas <- apply(fout$betaps, 2, sd)
  post.est$betam <- betam
  post.est$betas <- betas
  if (family == "negative.binomial" || family == "poisson.gamma") {
    kappam <- mean(fout$kappas)
    kappas <- sd(fout$kappas)
    post.est$kappam <- kappam
    post.est$kappas <- kappas
  }

  if (marginal.likelihood) {
    if (family == "negative.binomial" || family == "poisson.gamma") {
      betapost <- cbind(fout$betaps, fout$kappas)
    } else {
      betapost <- fout$betaps
    }
    beta_mn <- colMeans(betapost)
    beta_cov <- cov(betapost)
    beta_covi <- solve(beta_cov)
    lndetbcov <- log(det(beta_cov))

    logg <- .Fortran("gbglmgetlogg", as.matrix(betapost), as.integer(smcmc),
                     as.integer(ncol(betapost)), as.double(beta_mn),
                     as.matrix(beta_covi), as.double(lndetbcov),
                     logg = numeric(smcmc), NAOK = TRUE, PACKAGE = "bsamGP")$logg

    logjoint <- fout$loglikeps + fout$logpriorps
    ratiog <- logg - logjoint
    mratio <- max(ratiog)
    lilg <- exp(ratiog - mratio)
    lil <- -mratio - log(mean(lilg))
  }

  res.out <- list()
  res.out$call <- cl
  res.out$model <- "gblr"
  res.out$family <- family
  res.out$link <- link
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
  res.out$loglik.draws <- loglik.draws

  res.out$post.est <- post.est

  res.out$marglik <- marginal.likelihood
  if (marginal.likelihood) {
    res.out$lmarg <- lil[1]
  }

  if (family == 'bernoulli' && link == 'logit') {
    res.out$algorithm <- algorithm
  }

  class(res.out) <- "blm"
  res.out
}
