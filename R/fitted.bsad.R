"fitted.bsad" <- function(object, alpha = 0.05, HPD = TRUE, ...) {
  smcmc <- object$mcmc$smcmc

  fparg <- object$fit.draws$fpar
  fsemig <- object$fit.draws$fsemi
  fsemiMaxKappag <- object$fit.draws$fsemiMaxKappa

  if (object$parametric != "none") {
    fpar <- list()
    fparm <- apply(fparg, 2, mean)
    fpar$mean <- fparm
  }

  fsemi <- list()
  fsemim <- apply(fsemig, 2, mean)
  fsemi$mean <- fsemim

  fsemiMaxKappa <- list()
  fsemiMaxKappam <- apply(fsemiMaxKappag, 2, mean)
  fsemiMaxKappa$mean <- fsemiMaxKappam

  if (HPD) {
    n <- object$nint
    prob <- 1 - alpha

    if (object$parametric != "none") {
      fparg.o <- apply(fparg, 2, sort)
      gap <- max(1, min(smcmc - 1, round(smcmc * prob)))
      init <- 1:(smcmc - gap)
      inds <- apply(fparg.o[init + gap, , drop = FALSE] - fparg.o[init, , drop = FALSE], 2, which.min)
      fpar$lower <- fparg.o[cbind(inds, 1:n)]
      fpar$upper <- fparg.o[cbind(inds + gap, 1:n)]
    }

    fsemig.o <- apply(fsemig, 2, sort)
    gap <- max(1, min(smcmc - 1, round(smcmc * prob)))
    init <- 1:(smcmc - gap)
    inds <- apply(fsemig.o[init + gap, , drop = FALSE] - fsemig.o[init, , drop = FALSE], 2, which.min)
    fsemi$lower <- fsemig.o[cbind(inds, 1:n)]
    fsemi$upper <- fsemig.o[cbind(inds + gap, 1:n)]

    fsemiMaxKappag.o <- apply(fsemiMaxKappag, 2, sort)
    gap <- max(1, min(smcmc - 1, round(smcmc * prob)))
    init <- 1:(smcmc - gap)
    inds <- apply(fsemiMaxKappag.o[init + gap, , drop = FALSE] - fsemiMaxKappag.o[init, , drop = FALSE], 2, which.min)
    fsemiMaxKappa$lower <- fsemiMaxKappag.o[cbind(inds, 1:n)]
    fsemiMaxKappa$upper <- fsemiMaxKappag.o[cbind(inds + gap, 1:n)]
  } else {
    if (object$parametric != "none") {
      fpar$lower <- apply(fparg, 2, function(x) quantile(x, probs = alpha/2))
      fpar$upper <- apply(fparg, 2, function(x) quantile(x, probs = 1 - alpha/2))
    }

    fsemi$lower <- apply(fsemig, 2, function(x) quantile(x, probs = alpha/2))
    fsemi$upper <- apply(fsemig, 2, function(x) quantile(x, probs = 1 - alpha/2))

    fsemiMaxKappa$lower <- apply(fsemiMaxKappag, 2, function(x) quantile(x, probs = alpha/2))
    fsemiMaxKappa$upper <- apply(fsemiMaxKappag, 2, function(x) quantile(x, probs = 1 - alpha/2))
  }

  out <- object
  out$alpha <- alpha
  out$HPD <- HPD
  out$parametric <- object$parametric
  if (object$parametric != "none")
    out$fpar <- fpar
  out$fsemi <- fsemi
  out$fsemiMaxKappa <- fsemiMaxKappa
  class(out) <- "fitted.bsad"
  out
}
