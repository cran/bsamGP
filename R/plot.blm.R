"plot.blm" <- function(x, ...) {
  smcmc <- x$mcmc$smcmc
  param <- x$mcmc.draws$beta
  wname <- x$wname
  if (x$model != "gblr") {
    sigma2 <- x$mcmc.draws$sigma2
    param <- cbind(param, sigma2)
    wname <- c(wname, "sigma2")
  }
  p <- ncol(param)
  if (p == 2) {
    par(mfcol = c(2, 2))
    for (i in 1:p) {
      plot(1:smcmc, param[, i], xlab = "Iteration", ylab = "", main = wname[i], type = "l")
      plot(density(param[, i]), main = "")
    }
  } else if (p == 3) {
    par(mfcol = c(2, 3))
    for (i in 1:p) {
      plot(1:smcmc, param[, i], xlab = "Iteration", ylab = "", main = wname[i], type = "l")
      plot(density(param[, i]), main = "")
    }
  } else if (p == 4) {
    par(mfcol = c(2, 2))
    for (i in 1:p) {
      plot(1:smcmc, param[, i], xlab = "Iteration", ylab = "", main = wname[i], type = "l")
      plot(density(param[, i]), main = "")
      par(...)
    }
  } else {
    par(mfcol = c(2, 3))
    for (i in 1:p) {
      plot(1:smcmc, param[, i], xlab = "Iteration", ylab = "", main = wname[i], type = "l")
      plot(density(param[, i]), main = "")
      par(...)
    }
  }
}
