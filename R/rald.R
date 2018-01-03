"rald" <- function(n, location = 0, scale = 1, p = 0.5) {
  use.n <- if ((length.n <- length(n)) > 1)
    length.n else if (!is.numeric(n))
      stop("bad input for argument 'n'") else n
  tau <- p
  kappa <- sqrt(tau/(1 - tau))
  location <- rep(location, length.out = use.n)
  scale <- rep(scale, length.out = use.n)
  tau <- rep(tau, length.out = use.n)
  kappa <- rep(kappa, length.out = use.n)
  ans <- location + scale * log(runif(use.n)^kappa/runif(use.n)^(1/kappa))/sqrt(2)
  indexTF <- (scale > 0) & (tau > 0) & (tau < 1) & (kappa > 0)
  ans[!indexTF] <- NaN
  ans
}
