"predict.blm" <- function(object, newdata, alpha = 0.05, HPD = TRUE, ...) {
  smcmc <- object$mcmc$smcmc
  if (missing(newdata) || is.null(newdata)) {
    wbg <- object$fit.draws$wbeta
    newdata <- object$w
  } else {
    if (is.matrix(newdata)) {
      newdata <- as.matrix(newdata)
    }
    newdata <- cbind(1, newdata)
    colnames(newdata) <- object$wnames
    wbg <- object$mcmc.draws$beta %*% t(newdata)
  }
  n <- nrow(newdata)

  wbeta <- list()
  wbm <- apply(wbg, 2, mean)
  wbeta$mean <- wbm

  if (HPD) {
    prob <- 1 - alpha
    wbg.o <- apply(wbg, 2, sort)
    gap <- max(1, min(smcmc - 1, round(smcmc * prob)))
    init <- 1:(smcmc - gap)
    inds <- apply(wbg.o[init + gap, , drop = FALSE] - wbg.o[init, , drop = FALSE], 2, which.min)
    wbeta$lower <- wbg.o[cbind(inds, 1:n)]
    wbeta$upper <- wbg.o[cbind(inds + gap, 1:n)]
  } else {
    wbeta$lower <- apply(wbg, 2, function(x) quantile(x, prob = alpha/2))
    wbeta$upper <- apply(wbg, 2, function(x) quantile(x, prob = 1 - alpha/2))
  }

  if (object$model == 'gblr') {
    yhat <- list()
    if (object$link == 'log') {
      yhat$mean <- exp(wbeta$mean)
      yhat$upper <- exp(wbeta$upper)
      yhat$lower <- exp(wbeta$lower)
    } else if (object$link == 'probit') {
      yhat$mean <- pnorm(wbeta$mean)
      yhat$upper <- pnorm(wbeta$upper)
      yhat$lower <- pnorm(wbeta$lower)
    } else if (object$link == 'logit') {
      logit <- function(xx) 1 / (1 + exp(-xx))
      yhat$mean <- logit(wbeta$mean)
      yhat$upper <- logit(wbeta$upper)
      yhat$lower <- logit(wbeta$lower)
    }
  }

  out <- list()
  out$alpha <- alpha
  out$newdata <- newdata
  out$n <- n
  out$wbeta <- wbeta
  if (object$model == 'blq')
    out$p <- object$p
  if (object$model == 'gblr') {
    out$yhat <- yhat
    out$family <- object$family
    out$link <- object$link
  }
  class(out) <- "predict.blm"
  out
}
