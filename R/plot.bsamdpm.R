"plot.bsamdpm" <- function(x, ...) {
  nfun <- x$nfun
  nblow <- x$mcmc$nblow
  smcmc <- x$mcmc$smcmc
  betag <- x$mcmc.draws$beta
  sigmag <- x$mcmc.draws$sigma
  thetag <- x$mcmc.draws$theta
  taug <- x$mcmc.draws$tau
  gammag <- x$mcmc.draws$gamma
  zetag <- x$mcmc.draws$zeta

  t <- seq(nblow + 1, by = 1, length = smcmc)
  par(mfrow = c(3, 2))
  if (x$location) {
    mu0g <- x$dpm.draws$mu0
    plot(t, mu0g, type = "l", ylab = "", xlab = "", ylim = range(cbind(mu0g, betag)))
  } else {
    plot(t, betag[, 1], type = "l", ylab = "", xlab = "", ylim = range(betag))
  }
  if (ncol(betag) > 1) {
    for (i in 2:ncol(betag)) lines(t, betag[, i], col = i)
  }
  title("Beta")

  if (x$model != "bsardpm" && x$model != "bsaqdpm") {
    plot(t, sigmag, type = "l", ylab = "", xlab = "")
    title("Sigma")
  } else {
    plot(t, sigmag[, 1], type = "l", ylab = "", xlab = "", ylim = range(sigmag))
    for (i in 2:ncol(sigmag)) lines(t, sigmag[, i], col = i)
    title("Sigma")
  }

  for (i in 1:nfun) {
    plot(t, thetag[1, i, ], type = "l", ylab = "", xlab = "", ylim = range(thetag[, i, ]))
    for (j in 2:10) lines(t, thetag[j, i, ], col = j)
    title(paste("Function ", i, ": ", "Theta", sep = ""))

    plot(t, taug[, i], type = "l", ylab = "", xlab = "")
    title(paste("Function ", i, ": ", "Tau", sep = ""))

    plot(t, gammag[, i], type = "l", ylab = "", xlab = "")
    title(paste("Function ", i, ": ", "Gamma", sep = ""))

    plot(t, zetag[, i], type = "l", ylab = "", xlab = "")
    title(paste("Function ", i, ": ", "Zeta", sep = ""))

    if (nfun >= 2)
      par(mfrow = c(2, 2), ...)
  }

  par(mfrow = c(2, 2), ...)
  plot(t, x$dpm.draws$tmass, type = "l", ylab = "", xlab = "")
  title("Total Mass")

  plot(density(x$dpm.draws$tmass), lwd = 2, ylab = "", xlab = "", main = "")
  title("Total Mass")

  plot(t, x$dpm.draws$nclass, type = "l", ylab = "", xlab = "")
  title("Number of Clusters")

  barplot(table(x$dpm.draws$nclass)/smcmc, col = "gray", ylab = "", xlab = "")
  title("Number of Clusters")
}
