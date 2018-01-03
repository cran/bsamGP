"fitted.bsamdpm" <- function(object, alpha = 0.05, HPD = TRUE, ...) {
  smcmc <- object$mcmc$smcmc
  n <- object$n
  nint <- object$nint + 1
  nfun <- object$nfun

  edensg <- object$dpm.draws$edens
  fxobsg <- object$fit.draws$fxobs
  fxgridg <- object$fit.draws$fxgrid
  wbg <- object$fit.draws$wbeta
  yhatg <- object$fit.draws$yhat

  edens <- list()
  edensm <- apply(edensg, 2, mean)
  edens$mean <- edensm

  fxobs <- list()
  fxobsm <- apply(fxobsg, c(1, 2), mean)
  fxobs$mean <- fxobsm

  fxgrid <- list()
  fxgridm <- apply(fxgridg, c(1, 2), mean)
  fxgrid$mean <- fxgridm

  wbeta <- list()
  wbm <- apply(wbg, 2, mean)
  wbeta$mean <- wbm

  yhat <- list()
  ym <- apply(yhatg, 2, mean)
  yhat$mean <- ym

  if (HPD) {
    prob <- 1 - alpha

    edensg.o <- apply(edensg, 2, sort)
    gap <- max(1, min(smcmc - 1, round(smcmc * prob)))
    init <- 1:(smcmc - gap)
    inds <- apply(edensg.o[init + gap, , drop = FALSE] - edensg.o[init, , drop = FALSE], 2, which.min)
    edens$lower <- edensg.o[cbind(inds, 1:ncol(edensg))]
    edens$upper <- edensg.o[cbind(inds + gap, 1:ncol(edensg))]

    fx.l <- fx.u <- matrix(0, n, nfun)
    fxg.l <- fxg.u <- matrix(0, nint, nfun)
    for (i in 1:nfun) {
      fxobsg.o <- apply(fxobsg[, i, ], 1, sort)
      gap <- max(1, min(smcmc - 1, round(smcmc * prob)))
      init <- 1:(smcmc - gap)
      inds <- apply(fxobsg.o[init + gap, , drop = FALSE] - fxobsg.o[init, , drop = FALSE], 2, which.min)
      fx.l[, i] <- fxobsg.o[cbind(inds, 1:n)]
      fx.u[, i] <- fxobsg.o[cbind(inds + gap, 1:n)]

      fxgridg.o <- apply(fxgridg[, i, ], 1, sort)
      gap <- max(1, min(smcmc - 1, round(smcmc * prob)))
      init <- 1:(smcmc - gap)
      inds <- apply(fxgridg.o[init + gap, , drop = FALSE] - fxgridg.o[init, , drop = FALSE], 2, which.min)
      fxg.l[, i] <- fxgridg.o[cbind(inds, 1:nint)]
      fxg.u[, i] <- fxgridg.o[cbind(inds + gap, 1:nint)]
    }
    fxobs$lower <- fx.l
    fxobs$upper <- fx.u
    fxgrid$lower <- fxg.l
    fxgrid$upper <- fxg.u

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
    edens$lower <- apply(edensg, 2, function(x) quantile(x, prob = alpha/2))
    edens$upper <- apply(edensg, 2, function(x) quantile(x, prob = 1 - alpha/2))

    fxobs$lower <- apply(fxobsg, c(1, 2), function(x) quantile(x, prob = alpha/2))
    fxobs$upper <- apply(fxobsg, c(1, 2), function(x) quantile(x, prob = 1 - alpha/2))

    fxgrid$lower <- apply(fxgridg, c(1, 2), function(x) quantile(x, prob = alpha/2))
    fxgrid$upper <- apply(fxgridg, c(1, 2), function(x) quantile(x, prob = 1 - alpha/2))

    wbeta$lower <- apply(wbg, 2, function(x) quantile(x, prob = alpha/2))
    wbeta$upper <- apply(wbg, 2, function(x) quantile(x, prob = 1 - alpha/2))

    yhat$lower <- apply(yhatg, 2, function(x) quantile(x, prob = alpha/2))
    yhat$upper <- apply(yhatg, 2, function(x) quantile(x, prob = 1 - alpha/2))
  }

  out <- object
  out$alpha <- alpha
  out$HPD <- HPD
  out$yhat <- yhat
  out$wbeta <- wbeta
  out$fxobs <- fxobs
  out$fxgrid <- fxgrid
  out$xgrid <- object$fit.draws$xgrid
  out$edens <- edens
  out$egrid <- object$dpm.draws$egrid
  class(out) <- "fitted.bsamdpm"
  out
}
