"predict.bsamdpm" <- function(object, newp, newnp, alpha = 0.05, HPD = TRUE, ...) {
  smcmc <- object$mcmc$smcmc
  nbasis <- object$nbasis
  nint <- object$nint + 1
  nfun <- object$nfun
  fmodel <- object$fmodel
  fpm <- object$fpm
  xmin <- object$xmin
  xmax <- object$xmax

  if (missing(newp) && missing(newnp)) {
    n <- object$n
    newp <- object$w
    newnp <- object$x
    fxobsg <- object$fit.draws$fxobs
    wbg <- object$fit.draws$wbeta
    yhatg <- object$fit.draws$yhat
  } else if (missing(newp) && !missing(newnp)) {
    newp <- object$w
    if (!is.matrix(newnp))
      newnp <- as.matrix(newnp)
    n <- object$n
    if (n != nrow(newnp))
      stop('The number of observations for both parametric and nonparametric components must be same.')
    wbg <- object$fit.draws$wbeta
    fxobsg <- .Fortran("predictbsam", as.matrix(newnp), as.double(xmin), as.double(xmax),as.integer(n),
                       as.integer(nfun), as.integer(nbasis), as.integer(nint), as.integer(fmodel),
                       as.double(fpm), as.integer(smcmc), as.array(object$mcmc.draws$theta),
                       as.matrix(object$mcmc.draws$alpha), as.matrix(object$mcmc.draws$psi),
                       as.matrix(object$mcmc.draws$omega), fxobsg = array(0, dim = c(n, nfun, smcmc)),
                       NAOK = TRUE, PACKAGE = "bsamGP")$fxobsg
    yhatg <- wbg + t(apply(fxobsg, c(1,3), sum))
  } else if (!missing(newp) && missing(newnp)) {
    newnp <- object$x
    if (!is.matrix(newp))
      newp <- as.matrix(newp)
    newp <- cbind(1, newp)
    n <- object$n
    if (n != nrow(newp))
      stop('The number of observations for both parametric and nonparametric components must be same.')
    fxobsg <- object$fit.draws$fxobs
    wbg <- object$mcmc.draws$beta %*% t(newp)
    yhatg <- wbg + t(apply(fxobsg, c(1,3), sum))
  } else if (!missing(newp) && !missing(newnp)) {
    if (!is.matrix(newp))
      newp <- as.matrix(newp)
    newp <- cbind(1, newp)
    if (!is.matrix(newnp))
      newnp <- as.matrix(newnp)
    if (nrow(newp) != nrow(newnp))
      stop('The number of observations for both parametric and nonparametric components must be same.')
    n <- nrow(newp)
    wbg <- object$mcmc.draws$beta %*% t(newp)
    fxobsg <- .Fortran("predictbsam", as.matrix(newnp), as.double(xmin), as.double(xmax),as.integer(n),
                       as.integer(nfun), as.integer(nbasis), as.integer(nint), as.integer(fmodel),
                       as.double(fpm), as.integer(smcmc), as.array(object$mcmc.draws$theta),
                       as.matrix(object$mcmc.draws$alpha), as.matrix(object$mcmc.draws$psi),
                       as.matrix(object$mcmc.draws$omega), fxobsg = array(0, dim = c(n, nfun, smcmc)),
                       NAOK = TRUE, PACKAGE = "bsamGP")$fxobsg
    yhatg <- wbg + t(apply(fxobsg, c(1,3), sum))
  }

  fxobs <- list()
  fxobsm <- apply(fxobsg, c(1, 2), mean)
  fxobs$mean <- fxobsm

  wbeta <- list()
  wbm <- apply(wbg, 2, mean)
  wbeta$mean <- wbm

  yhat <- list()
  ym <- apply(yhatg, 2, mean)
  yhat$mean <- ym

  if (HPD) {
    prob <- 1 - alpha

    fx.l <- fx.u <- matrix(0, n, nfun)
    for (i in 1:nfun) {
      fxobsg.o <- apply(fxobsg[, i, ], 1, sort)
      gap <- max(1, min(smcmc - 1, round(smcmc * prob)))
      init <- 1:(smcmc - gap)
      inds <- apply(fxobsg.o[init + gap, , drop = FALSE] - fxobsg.o[init, , drop = FALSE], 2, which.min)
      fx.l[, i] <- fxobsg.o[cbind(inds, 1:n)]
      fx.u[, i] <- fxobsg.o[cbind(inds + gap, 1:n)]
    }
    fxobs$lower <- fx.l
    fxobs$upper <- fx.u

    wbg.o <- apply(wbg, 2, sort)
    gap <- max(1, min(smcmc - 1, round(smcmc * prob)))
    init <- 1:(smcmc - gap)
    inds <- apply(wbg.o[init + gap, , drop = FALSE] - wbg.o[init, , drop = FALSE], 2, which.min)
    wbeta$lower <- wbg.o[cbind(inds, 1:n)]
    wbeta$upper <- wbg.o[cbind(inds + gap, 1:n)]

    yhatg.o <- apply(yhatg, 2, sort)
    gap <- max(1, min(smcmc - 1, round(smcmc * prob)))
    init <- 1:(smcmc - gap)
    inds <- apply(yhatg.o[init + gap, , drop = FALSE] - yhatg.o[init, , drop = FALSE], 2, which.min)
    yhat$lower <- yhatg.o[cbind(inds, 1:n)]
    yhat$upper <- yhatg.o[cbind(inds + gap, 1:n)]
  } else {
    fxobs$lower <- apply(fxobsg, c(1, 2), function(x) quantile(x, prob = alpha/2))
    fxobs$upper <- apply(fxobsg, c(1, 2), function(x) quantile(x, prob = 1 - alpha/2))

    wbeta$lower <- apply(wbg, 2, function(x) quantile(x, prob = alpha/2))
    wbeta$upper <- apply(wbg, 2, function(x) quantile(x, prob = 1 - alpha/2))

    yhat$lower <- apply(yhatg, 2, function(x) quantile(x, prob = alpha/2))
    yhat$upper <- apply(yhatg, 2, function(x) quantile(x, prob = 1 - alpha/2))
  }

  out <- list()
  out$n <- n
  out$nbasis <- nbasis
  out$newp <- newp
  out$newnp <- newnp
  out$alpha <- alpha
  out$HPD <- HPD
  out$yhat <- yhat
  out$wbeta <- wbeta
  out$fxobs <- fxobs
  if (object$model == 'bsaqdpm')
    out$p <- out$p
  class(out) <- "predict.bsamdpm"
  out
}
