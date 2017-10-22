"predict.bsam" <- function(object, newp, newnp, alpha = 0.05, HPD = TRUE, type = "response", ...) {
  if(type!="response" && type !="mean")
    stop("type has to be either 'response' or 'mean'")
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
    if (object$model != "gbsar") {
      yhatg <- object$fit.draws$yhat
    } else {
      yhatg <- object$fit.draws$muhat
    }
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
    if (object$model == "gbsar") {
      if (object$link == 'probit') {
        yhatg <- pnorm(wbg + t(apply(fxobsg, c(1,3), sum)))
      } else if (object$link == 'logit') {
        logit <- function(xx) 1 / (1 + exp(-xx))
        yhatg <- logit(wbg + t(apply(fxobsg, c(1,3), sum)))
      } else {
        yhatg <- exp(wbg + t(apply(fxobsg, c(1,3), sum)))
      }
    } else if (object$model == "bsaq") {
    	if(type=="response") {
    		p <- object$p
    		eta1 <- (1 - 2*p) / (p * (1-p))
    		eta22 <- 2/(p*(1-p))
    		zg <- t(sapply(object$mcmc.draws$sigma, function(sigma) rexp(n, 1/sigma)))
    		yhatg <- wbg + t(apply(fxobsg, c(1,3), sum)) + rnorm(n*smcmc, mean=eta1*zg, sd=sqrt(eta22*zg*object$mcmc.draws$sigma)) # Add error
    	} else {	# mean
    		yhatg <- wbg + t(apply(fxobsg, c(1,3), sum))
    	}
    } else {	# bsar
    	if(type=="response") {
    		yhatg <- wbg + t(apply(fxobsg, c(1,3), sum)) + t(sapply(object$mcmc.draws$sigma, function(sigma) rnorm(n, sd=sigma))) # Add error
    	} else {	# mean
    		yhatg <- wbg + t(apply(fxobsg, c(1,3), sum))
    	}
    }
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

    if (object$model == "gbsar") {
      if (object$link == 'probit') {
        yhatg <- pnorm(wbg + t(apply(fxobsg, c(1,3), sum)))
      } else if (object$link == 'logit') {
        logit <- function(xx) 1 / (1 + exp(-xx))
        yhatg <- logit(wbg + t(apply(fxobsg, c(1,3), sum)))
      } else {
        yhatg <- exp(wbg + t(apply(fxobsg, c(1,3), sum)))
      }
    } else if (object$model == "bsaq") {
    	if(type=="response") {
    		p <- object$p
    		eta1 <- (1 - 2*p) / (p * (1-p))
    		eta22 <- 2/(p*(1-p))
    		zg <- t(sapply(object$mcmc.draws$sigma, function(sigma) rexp(n, 1/sigma)))
    		yhatg <- wbg + t(apply(fxobsg, c(1,3), sum)) + rnorm(n*smcmc, mean=eta1*zg, sd=sqrt(eta22*zg*object$mcmc.draws$sigma)) # Add error
    	} else {	# mean
    		yhatg <- wbg + t(apply(fxobsg, c(1,3), sum))
      	}
    } else {	# bsar
    	if(type=="response") {
    		yhatg <- wbg + t(apply(fxobsg, c(1,3), sum)) + t(sapply(object$mcmc.draws$sigma, function(sigma) rnorm(n, sd=sigma))) # Add error
    	} else {	# mean
    		yhatg <- wbg + t(apply(fxobsg, c(1,3), sum))
    	}
    }
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

    if (object$model == "gbsar") {
      	if (object$link == 'probit') {
      	  yhatg <- pnorm(wbg + t(apply(fxobsg, c(1,3), sum)))
      	} else if (object$link == 'logit') {
      	  logit <- function(xx) 1 / (1 + exp(-xx))
      	  yhatg <- logit(wbg + t(apply(fxobsg, c(1,3), sum)))
      	} else {
      	  yhatg <- exp(wbg + t(apply(fxobsg, c(1,3), sum)))
      	}
    } else if (object$model == "bsaq") {
    	if(type=="response") {
    		p <- object$p
    		eta1 <- (1 - 2*p) / (p * (1-p))
    		eta22 <- 2/(p*(1-p))
    		zg <- t(sapply(object$mcmc.draws$sigma, function(sigma) rexp(n, 1/sigma)))
    		yhatg <- wbg + t(apply(fxobsg, c(1,3), sum)) + rnorm(n*smcmc, mean=eta1*zg, sd=sqrt(eta22*zg*object$mcmc.draws$sigma)) # Add error
    	} else {	# mean
    		yhatg <- wbg + t(apply(fxobsg, c(1,3), sum))
      }
    } else {	# bsar
    	if(type=="response") {
    		yhatg <- wbg + t(apply(fxobsg, c(1,3), sum)) + t(sapply(object$mcmc.draws$sigma, function(sigma) rnorm(n, sd=sigma))) # Add error
    	} else {	# mean
    		yhatg <- wbg + t(apply(fxobsg, c(1,3), sum))
    	}
    }
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

  ## Parametric Residual ##
  fxResidg <- array(0, c(n,nfun,smcmc))
  for (i in 1:nfun) {
    fxResidg[,i,] <- t(yhatg - wbg) - apply(fxobsg[,-i,,drop=FALSE], c(1,3), sum)
  }
  
  fxResid <- list()
  fxResid$mean <- apply(fxResidg, c(1, 2), mean)

  if (HPD) {
    prob <- 1 - alpha

    fx.l <- fx.u <- matrix(0, n, nfun)
    gap <- max(1, min(smcmc - 1, round(smcmc * prob)))
    init <- 1:(smcmc - gap)
    for (i in 1:nfun) {
      fxobsg.o <- apply(fxobsg[, i, ], 1, sort)
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

    ## Parametric Residual ##
    fxResid.l <- fxResid.u <- matrix(0, n, nfun)
    gap <- max(1, min(smcmc - 1, round(smcmc * prob)))
    init <- 1:(smcmc - gap)
    for (i in 1:nfun) {
      fxResidg.o <- apply(fxResidg[, i, ], 1, sort)
      inds <- apply(fxResidg.o[init + gap, , drop = FALSE] - fxResidg.o[init, , drop = FALSE], 2, which.min)
      fxResid.l[, i] <- fxResidg.o[cbind(inds, 1:n)]
      fxResid.u[, i] <- fxResidg.o[cbind(inds + gap, 1:n)]
    }
    fxResid$lower <- fxResid.l
    fxResid$upper <- fxResid.u
    ##
  } else {
    fxobs$lower <- apply(fxobsg, c(1, 2), function(x) quantile(x, prob = alpha/2))
    fxobs$upper <- apply(fxobsg, c(1, 2), function(x) quantile(x, prob = 1 - alpha/2))

    wbeta$lower <- apply(wbg, 2, function(x) quantile(x, prob = alpha/2))
    wbeta$upper <- apply(wbg, 2, function(x) quantile(x, prob = 1 - alpha/2))

    yhat$lower <- apply(yhatg, 2, function(x) quantile(x, prob = alpha/2))
    yhat$upper <- apply(yhatg, 2, function(x) quantile(x, prob = 1 - alpha/2))

    ## Parametric Residual ##
    fxResid$lower <- apply(fxResidg, c(1, 2), function(x) quantile(x, prob = alpha/2))
    fxResid$upper <- apply(fxResidg, c(1, 2), function(x) quantile(x, prob = 1 - alpha/2))
    ##    
  }

  out <- list()
  out$n <- n
  out$nbasis <- nbasis
  out$newp <- newp
  out$newnp <- newnp
  out$alpha <- alpha
  out$HPD <- HPD
  out$type <- type
  out$yhat <- yhat
  out$wbeta <- wbeta
  out$fxobs <- fxobs
  if (object$model == 'bsaq')
    out$p <- out$p
  if (object$model == 'gbsar') {
    out$family <- object$family
    out$link <- object$link
  } else {
    out$fxResid <- fxResid
  }

  class(out) <- "predict.bsam"
  out
}
