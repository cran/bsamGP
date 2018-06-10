"plot.fitted.bsamdpm" <- function(x, ask = FALSE, ggplot2 = TRUE, legend.position = "none", ...) {
  yobs <- x$y
  xobs <- x$x
  nobs <- x$n
  nfun <- x$nfun
  smcmc <- x$mcmc$smcmc
  xgrid <- x$fit.draws$xgrid
  ngrid <- x$nint + 1
  if (x$model == "bsardpm" && x$location) {
    mum <- x$post.est$mum
  }
  wbm <- x$wbeta$mean
  fxobsm <- x$fxobs$mean
  fxobsl <- x$fxobs$lower
  fxobsu <- x$fxobs$upper
  fxgridm <- x$fxgrid$mean
  fxgridl <- x$fxgrid$lower
  fxgridu <- x$fxgrid$upper
  prob <- (1 - x$alpha) * 100
  HPD <- x$HPD

  xname <- x$xname
  shape <- x$shape

  par(ask = ask, ...)
  if (ggplot2) {
    if (nfun == 1) {
      resid <- yobs - wbm
      if (x$model == "bsardpm" && x$location)
        resid <- resid - mum
      datp <- data.frame(x = xobs[, 1], y = resid, Estimates = rep("Parametric Residuals", nobs))
      if (HPD) {
        datl <- data.frame(x = rep(xgrid, 3), fx = c(fxgridu, fxgridm, fxgridl),
                           Estimates = c(rep(paste(prob, "% HPD UCI", sep = ""), ngrid),
                                         rep("Posterior Mean", ngrid),
                                         rep(paste(prob, "% HPD LCI", sep = ""), ngrid)))

        dato <- data.frame(x = rep(xobs, 3), fx = c(fxobsu, fxobsm, fxobsl),
                           Estimates = c(rep(paste(prob, "% HPD UCI", sep = ""), nobs),
                                         rep("Posterior Mean", nobs),
                                         rep(paste(prob, "% HPD LCI", sep = ""), nobs)))
      } else {
        datl <- data.frame(x = rep(xgrid, 3), fx = c(fxgridu, fxgridm, fxgridl),
                           Estimates = c(rep(paste(prob, "% equal-tail UCI", sep = ""), ngrid),
                                         rep("Posterior Mean", ngrid),
                                         rep(paste(prob, "% equal-tail LCI", sep = ""), ngrid)))

        dato <- data.frame(x = rep(xobs, 3), fx = c(fxobsu, fxobsm, fxobsl),
                           Estimates = c(rep(paste(prob, "% equal-tail UCI", sep = ""), nobs),
                                         rep("Posterior Mean", nobs),
                                         rep(paste(prob, "% equal-tail LCI", sep = ""), nobs)))
      }
      plt1 <- ggplot(datl)
      plt1 <- plt1 + aes_string(x = 'x', y = 'fx')
      plt1 <- plt1 + aes_string(group = 'Estimates')
      plt1 <- plt1 + aes_string(shape = 'Estimates', linetype = 'Estimates', colour = 'Estimates')
      plt1 <- plt1 + geom_line(size = 0.8)
      plt1 <- plt1 + xlab(x$xname[1])
      plt1 <- plt1 + ylab(paste('f(', x$xname[1], ')', sep = ''))
      plt1 <- plt1 + theme_bw()
      plt1 <- plt1 + theme(legend.position = legend.position)
      plt1 <- plt1 + scale_linetype_manual(values = c("dotdash", "dotdash", "solid"))

      plt2 <- ggplot(dato)
      plt2 <- plt2 + geom_point(data = datp, mapping = aes_string(x = 'x', y = 'y'), shape = 21, alpha = 0.6)
      plt2 <- plt2 + aes_string(x = 'x', y = 'fx')
      plt2 <- plt2 + aes_string(group = 'Estimates')
      plt2 <- plt2 + aes_string(shape = 'Estimates', linetype = 'Estimates', colour = 'Estimates')
      plt2 <- plt2 + geom_line(size = 0.8)
      plt2 <- plt2 + xlab(x$xname[1])
      if (x$model == "bsardpm" && x$location) {
        plt2 <- plt2 + ylab(expression(y - w^T * hat(beta) - hat(mu)))
      } else {
        plt2 <- plt2 + ylab(expression(y - w^T * hat(beta)))
      }
      plt2 <- plt2 + theme_bw()
      plt2 <- plt2 + theme(legend.position = legend.position)
      plt2 <- plt2 + scale_linetype_manual(values = c("dotdash", "dotdash", "dotted", "solid"))

      grid.arrange(plt1, plt2, nrow = 2)
    } else {
      resid <- yobs - wbm
      if (x$model == "bsardpm" && x$location)
        resid <- resid - mum
      for (i in 1:nfun) {
        datp <- data.frame(x = xobs[, i], y = resid - rowSums(fxobsm[, -i, drop = FALSE]),
                           Estimates = rep("Partial Residuals", nobs))
        if (HPD) {
          datl <- data.frame(x = rep(xgrid[, i], 3), fx = c(fxgridu[, i], fxgridm[, i], fxgridl[, i]),
                             Estimates = c(rep(paste(prob, "% HPD UCI", sep = ""), ngrid),
                                           rep("Posterior Mean", ngrid),
                                           rep(paste(prob, "% HPD LCI", sep = ""), ngrid)))

          dato <- data.frame(x = rep(xobs[, i], 3), fx = c(fxobsu[, i], fxobsm[, i], fxobsl[, i]),
                             Estimates = c(rep(paste(prob, "% HPD UCI", sep = ""), nobs),
                                           rep("Posterior Mean", nobs),
                                           rep(paste(prob, "% HPD LCI", sep = ""), nobs)))
        } else {
          datl <- data.frame(x = rep(xgrid[, i], 3), fx = c(fxgridu[, i], fxgridm[, i], fxgridl[, i]),
                             Estimates = c(rep(paste(prob, "% equal-tail UCI", sep = ""), ngrid),
                                           rep("Posterior Mean", ngrid),
                                           rep(paste(prob, "% equal-tail LCI", sep = ""), ngrid)))

          dato <- data.frame(x = rep(xobs[, i], 3), fx = c(fxobsu[, i], fxobsm[, i], fxobsl[, i]),
                             Estimates = c(rep(paste(prob, "% equal-tail UCI", sep = ""), nobs),
                                           rep("Posterior Mean", nobs),
                                           rep(paste(prob, "% equal-tail LCI", sep = ""), nobs)))
        }
        plt1 <- ggplot(datl)
        plt1 <- plt1 + aes_string(x = 'x', y = 'fx')
        plt1 <- plt1 + aes_string(group = 'Estimates')
        plt1 <- plt1 + aes_string(shape = 'Estimates', linetype = 'Estimates', colour = 'Estimates')
        plt1 <- plt1 + geom_line(size = 0.8)
        plt1 <- plt1 + xlab(parse(text=x$xname[i]))
        plt1 <- plt1 + ylab(paste('f(', x$xname[i], ')', sep = ''))
        plt1 <- plt1 + theme_bw()
        plt1 <- plt1 + theme(legend.position = legend.position)
        plt1 <- plt1 + scale_linetype_manual(values = c("dotdash", "dotdash", "solid"))

        plt2 <- ggplot(dato)
        plt2 <- plt2 + geom_point(data = datp, mapping = aes_string(x = 'x', y = 'y'), shape = 21, alpha = 0.6)
        plt2 <- plt2 + aes_string(x = 'x', y = 'fx')
        plt2 <- plt2 + aes_string(group = 'Estimates')
        plt2 <- plt2 + aes_string(shape = 'Estimates', linetype = 'Estimates', colour = 'Estimates')
        plt2 <- plt2 + geom_line(size = 0.8)
        plt2 <- plt2 + xlab(parse(text=x$xname[i]))
        plt2 <- plt2 + ylab("Partial Residuals")
        plt2 <- plt2 + theme_bw()
        plt2 <- plt2 + theme(legend.position = legend.position)
        plt2 <- plt2 + scale_linetype_manual(values = c("dotdash", "dotdash", "dotted", "solid"))

        grid.arrange(plt1, plt2, nrow = 2)
      }
    }

    edensm <- x$edens$mean
    edensl <- x$edens$lower
    edensu <- x$edens$upper
    egrid <- x$egrid
    negrid <- x$dpm.draws$ngrid
    if (HPD) {
      datl <- data.frame(x = rep(egrid, 3), fx = c(edensu, edensm, edensl),
                         Estimates = c(rep(paste(prob, "% HPD UCI", sep = ""), negrid),
                                       rep("Posterior Mean", negrid),
                                       rep(paste(prob, "% HPD LCI", sep = ""), negrid)))
    } else {
      datl <- data.frame(x = rep(egrid, 3), fx = c(edensu, edensm, edensl),
                         Estimates = c(rep(paste(prob, "% equal-tail UCI", sep = ""), negrid),
                                       rep("Posterior Mean", negrid),
                                       rep(paste(prob, "% equal-tail LCI", sep = ""), negrid)))
    }
    plt3 <- ggplot(datl)
    plt3 <- plt3 + aes_string(x = 'x', y = 'fx')
    plt3 <- plt3 + aes_string(group = 'Estimates')
    plt3 <- plt3 + aes_string(shape = 'Estimates', linetype = 'Estimates', colour = 'Estimates')
    plt3 <- plt3 + geom_line(size = 0.8)
    plt3 <- plt3 + xlab(expression(epsilon))
    plt3 <- plt3 + ylab("Density")
    plt3 <- plt3 + theme_bw()
    plt3 <- plt3 + theme(legend.position = legend.position)
    plt3 <- plt3 + scale_linetype_manual(values = c("dotdash", "dotdash", "solid"))
    grid.arrange(plt3, nrow = 1)
  } else {
    resid <- yobs - wbm
    if (x$model == "bsardpm" && x$location)
      resid <- resid - mum
    if (nfun == 1) {
      plot(x = xgrid, y = fxgridm, pch = NA, main = '',
           ylim = range(c(fxgridl, fxgridu)),
           xlab = x$xname[1], ylab = paste('f(', x$xname[1], ')', sep = ''), ...)
      polygon(x = c(xgrid, rev(xgrid)),
              y = c(fxgridl, rev(fxgridu)), col = 'gray70', lty = 2)
      lines(x = xgrid, y = fxgridm, lwd = 2, col = 2)

      plot(x = xobs, y = resid, pch = NA, main = '',
           xlab = x$xname[1], ylab = 'Parametric Residuals', ...)
      polygon(x = c(xobs, rev(xobs)),
              y = c(fxobsl, rev(fxobsu)), col = 'gray70', lty = 2)
      points(x = xobs, y = resid, lwd = 2)
      lines(x = xobs, y = fxobsm, lwd = 2, col = 2)
    } else {
      for (i in 1:nfun) {
        plot(x = xgrid[, i], y = fxgridm[, i], pch = NA,
             ylim = range(c(fxgridl[, i], fxgridu[, i])), main = '',
             xlab = x$xname[i], ylab = paste('f(', x$xname[i], ')', sep = ''), ...)
        polygon(x = c(xgrid[, i], rev(xgrid[, i])),
                y = c(fxgridl[, i], rev(fxgridu[, i])), col = 'gray70', lty = 2)
        lines(x = xgrid[, i], y = fxgridm[, i], lwd = 2, col = 2)

        plot(x = xobs[, i], y = resid - rowSums(fxobsm[, -i, drop = FALSE]), pch = NA,
             main = '', xlab = x$xname[i], ylab = 'Partial Residuals', ...)
        polygon(x = c(xobs[, i], rev(xobs[, i])),
                y = c(fxobsl[, i], rev(fxobsu[, i])), col = 'gray70', lty = 2)
        points(x = xobs[, i], y = resid - rowSums(fxobsm[, -i, drop = FALSE]), lwd = 2)
        lines(x = xobs[, i], y = fxobsm[, i], lwd = 2, col = 2)
      }
    }

    edensm <- x$edens$mean
    edensl <- x$edens$lower
    edensu <- x$edens$upper
    egrid <- x$egrid
    negrid <- x$dpm.draws$ngrid
    plot(egrid, edensm, type = 'l', lwd = 2, lty = 1, ylim = c(0, max(edensu)),
         xlab = expression(epsilon), ylab = 'Density', col = "dodgerblue", ...)
    polygon(x = c(egrid, rev(egrid)), y = c(edensl, rev(edensu)),
            col = 'gray70', lty = 2)
    lines(egrid, edensm, lwd = 2, col = 2)
  }
}
