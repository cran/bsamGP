"plot.bsad" <- function(x, ...) {
  nblow <- x$mcmc$nblow
  smcmc <- x$mcmc$smcmc
  betag <- x$mcmc.draws$beta
  sigmag <- x$mcmc.draws$sigma
  thetag <- x$mcmc.draws$theta
  taug <- x$mcmc.draws$tau
  gammag <- x$mcmc.draws$gam
  kappag <- x$mcmc.draws$kappa

  t <- seq(nblow + 1, by = 1, length = smcmc)
  par(mfrow = c(3, 2))
  plot(t, betag[, 1], type = "l", ylab = "", xlab = "", ylim = range(betag))
  if (ncol(betag) > 1) {
    for (i in 2:ncol(betag)) lines(t, betag[, i], col = i)
  }
  title("Beta")

  plot(t, sigmag, type = "l", ylab = "", xlab = "")
  title("Sigma")

  plot(t, thetag[, 1], type = "l", ylab = "", xlab = "", ylim = range(thetag))
  for (j in 2:10) lines(t, thetag[, j], col = j)
  title("Theta")

  plot(t, taug, type = "l", ylab = "", xlab = "")
  title("Tau")

  plot(t, gammag, type = "l", ylab = "", xlab = "")
  title("Gamma")

  plot(t, kappag, type = "l", ylab = "", xlab = "")
  title("Kappa")
}
