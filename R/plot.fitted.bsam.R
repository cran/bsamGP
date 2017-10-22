"plot.fitted.bsam" <- function(x, ask = FALSE, ggplot2 = TRUE, ...) {
  yobs <- x$y
  xobs <- x$x
  nobs <- x$n
  nfun <- x$nfun
  smcmc <- x$mcmc$smcmc
  xgrid <- x$fit.draws$xgrid
  ngrid <- x$nint + 1
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
                           Estimates = c(rep(paste(prob, "% Equal-tail UCI", sep = ""), ngrid),
                                         rep("Posterior Mean", ngrid),
                                         rep(paste(prob, "% Equal-tail LCI", sep = ""), ngrid)))
        
        dato <- data.frame(x = rep(xobs, 3), fx = c(fxobsu, fxobsm, fxobsl),
                           Estimates = c(rep(paste(prob, "% Equal-tail UCI", sep = ""), nobs),
                                         rep("Posterior Mean", nobs),
                                         rep(paste(prob, "% Equal-tail LCI", sep = ""), nobs)))
      }
      plt1 <- ggplot(datl)
      plt1 <- plt1 + aes_string(x = 'x', y = 'fx')
      plt1 <- plt1 + aes_string(group = 'Estimates')
      plt1 <- plt1 + aes_string(shape = 'Estimates', linetype = 'Estimates', colour = 'Estimates')
      plt1 <- plt1 + geom_line(size = 0.8)
      plt1 <- plt1 + xlab(x$xname[1])
      plt1 <- plt1 + ylab(paste('f(', x$xname[1], ')', sep = ''))
      plt1 <- plt1 + theme_bw()
      plt1 <- plt1 + theme(legend.position = "top")
      plt1 <- plt1 + scale_linetype_manual(values = c("dotdash", "dotdash", "solid"))
      
      if (x$model == "gbsar") {
        if (x$family == "bernoulli") {
          if (x$link == 'probit') {
            dato$fx = pnorm(c(rep(median(x$wbeta$upper), nobs), 
                              rep(median(x$wbeta$mean), nobs), 
                              rep(median(x$wbeta$lower), nobs)) + dato$fx)
          } else {
            logit <- function(xx) 1 / (1 + exp(-xx))
            dato$fx = logit(c(rep(median(x$wbeta$upper), nobs), 
                              rep(median(x$wbeta$mean), nobs), 
                              rep(median(x$wbeta$lower), nobs)) + dato$fx)
          }
          datp <- data.frame(x = xobs[, 1], y = yobs, Estimates = rep("Observations", nobs))
          plt2 <- ggplot(dato)
          plt2 <- plt2 + geom_point(data = datp, mapping = aes_string(x = 'x', y = 'y'), shape = 21, alpha = 0.6)
          plt2 <- plt2 + aes_string(x = 'x', y = 'fx')
          plt2 <- plt2 + aes_string(group = 'Estimates')
          plt2 <- plt2 + aes_string(shape = 'Estimates', linetype = 'Estimates', colour = 'Estimates')
          plt2 <- plt2 + geom_line(size = 0.8)
          plt2 <- plt2 + xlab(x$xname[1])
          plt2 <- plt2 + theme_bw()
          plt2 <- plt2 + theme(legend.position = "top")
          plt2 <- plt2 + ylab(paste('P(', x$yname[1], ')', sep = ''))
          plt2 <- plt2 + ylim(c(0,1))
          plt2 <- plt2 + scale_linetype_manual(values = c("dotdash", "dotdash", "solid", "solid"))
        } else {
          dato$fx = exp(c(rep(median(x$wbeta$upper), nobs), 
                          rep(median(x$wbeta$mean), nobs), 
                          rep(median(x$wbeta$lower), nobs)) + dato$fx)
          datp <- data.frame(x = xobs[, 1], y = yobs, Estimates = rep("Observations", nobs))
          plt2 <- ggplot(dato)
          plt2 <- plt2 + geom_point(data = datp, mapping = aes_string(x = 'x', y = 'y'), shape = 21, alpha = 0.6)
          plt2 <- plt2 + aes_string(x = 'x', y = 'fx')
          plt2 <- plt2 + aes_string(group = 'Estimates')
          plt2 <- plt2 + aes_string(shape = 'Estimates', linetype = 'Estimates', colour = 'Estimates')
          plt2 <- plt2 + geom_line(size = 0.8)
          plt2 <- plt2 + xlab(x$xname[1])
          plt2 <- plt2 + theme_bw()
          plt2 <- plt2 + theme(legend.position = "top")
          plt2 <- plt2 + ylab(x$yname[1])
          plt2 <- plt2 + scale_linetype_manual(values = c("dotdash", "dotdash", "solid", "solid"))
        }
      } else {
        datp <- data.frame(x = xobs[, 1], y = yobs - wbm, Estimates = rep("Parametric Residuals", nobs))
        plt2 <- ggplot(dato)
        plt2 <- plt2 + geom_point(data = datp, mapping = aes_string(x = 'x', y = 'y'), shape = 21, alpha = 0.6)
        plt2 <- plt2 + aes_string(x = 'x', y = 'fx')
        plt2 <- plt2 + aes_string(group = 'Estimates')
        plt2 <- plt2 + aes_string(shape = 'Estimates', linetype = 'Estimates', colour = 'Estimates')
        plt2 <- plt2 + geom_line(size = 0.8)
        plt2 <- plt2 + xlab(x$xname[1])
        plt2 <- plt2 + theme_bw()
        plt2 <- plt2 + theme(legend.position = "top")
        plt2 <- plt2 + ylab('Parametric Residuals')
        plt2 <- plt2 + scale_linetype_manual(values = c("dotdash", "dotdash", "solid", "solid"))
      }
      grid.arrange(plt1, plt2, nrow = 2)
    } else {
      for (i in 1:nfun) {
        if (!ask)
          dev.new()
        
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
                             Estimates = c(rep(paste(prob, "% Equal-tail UCI", sep = ""), ngrid),
                                           rep("Posterior Mean", ngrid),
                                           rep(paste(prob, "% Equal-tail LCI", sep = ""), ngrid)))
          
          dato <- data.frame(x = rep(xobs[, i], 3), fx = c(fxobsu[, i], fxobsm[, i], fxobsl[, i]),
                             Estimates = c(rep(paste(prob, "% Equal-tail UCI", sep = ""), nobs),
                                           rep("Posterior Mean", nobs),
                                           rep(paste(prob, "% Equal-tail LCI", sep = ""), nobs)))
        }
        plt1 <- ggplot(datl)
        plt1 <- plt1 + aes_string(x = 'x', y = 'fx')
        plt1 <- plt1 + aes_string(group = 'Estimates')
        plt1 <- plt1 + aes_string(shape = 'Estimates', linetype = 'Estimates', colour = 'Estimates')
        plt1 <- plt1 + geom_line(size = 0.8)
        plt1 <- plt1 + xlab(parse(text=x$xname[i]))
        plt1 <- plt1 + ylab(paste('f(', x$xname[i], ')', sep = ''))
        plt1 <- plt1 + theme_bw()
        plt1 <- plt1 + theme(legend.position = "top")
        plt1 <- plt1 + scale_linetype_manual(values = c("dotdash", "dotdash", "solid"))
        
        if (x$model == "gbsar") {
          if (x$family == "bernoulli") {
            if (x$link == "probit") {
              dato$fx = pnorm(dato$fx)
            } else {
              logit <- function(xx) 1 / (1 + exp(-xx))
              dato$fx = logit(dato$fx)
            }
            datp <- data.frame(x = xobs[, i], y = yobs, Estimates = rep("Observations", nobs))
            plt2 <- ggplot(dato)
            plt2 <- plt2 + geom_point(data = datp, mapping = aes_string(x = 'x', y = 'y'), shape = 21, alpha = 0.6)
            plt2 <- plt2 + aes_string(x = 'x', y = 'fx')
            plt2 <- plt2 + aes_string(group = 'Estimates')
            plt2 <- plt2 + aes_string(shape = 'Estimates', linetype = 'Estimates', colour = 'Estimates')
            plt2 <- plt2 + geom_line(size = 0.8)
            plt2 <- plt2 + xlab(x$xname[i])
            plt2 <- plt2 + theme_bw()
            plt2 <- plt2 + theme(legend.position = "top")
            plt2 <- plt2 + ylab(paste('P(', x$yname[1], ')', sep = ''))
            plt2 <- plt2 + ylim(c(0,1))
            plt2 <- plt2 + scale_linetype_manual(values = c("dotdash", "dotdash", "solid", "solid"))
          } else {
            dato$fx = exp(dato$fx)
            datp <- data.frame(x = xobs[, i], y = yobs, Estimates = rep("Observations", nobs))
            plt2 <- ggplot(dato)
            plt2 <- plt2 + geom_point(data = datp, mapping = aes_string(x = 'x', y = 'y'), shape = 21, alpha = 0.6)
            plt2 <- plt2 + aes_string(x = 'x', y = 'fx')
            plt2 <- plt2 + aes_string(group = 'Estimates')
            plt2 <- plt2 + aes_string(shape = 'Estimates', linetype = 'Estimates', colour = 'Estimates')
            plt2 <- plt2 + geom_line(size = 0.8)
            plt2 <- plt2 + xlab(x$xname[i])
            plt2 <- plt2 + theme_bw()
            plt2 <- plt2 + theme(legend.position = "top")
            plt2 <- plt2 + ylab(x$yname[1])
            plt2 <- plt2 + scale_linetype_manual(values = c("dotdash", "dotdash", "solid", "solid"))
          }
        } else {
          datp <- data.frame(x = xobs[, i], y = yobs - wbm - rowSums(fxobsm[, -i, drop = FALSE]),
                             Estimates = rep("Partial Residuals", nobs))
          
          plt2 <- ggplot(dato)
          plt2 <- plt2 + geom_point(data = datp, mapping = aes_string(x = 'x', y = 'y'), shape = 21, alpha = 0.6)
          plt2 <- plt2 + aes_string(x = 'x', y = 'fx')
          plt2 <- plt2 + aes_string(group = 'Estimates')
          plt2 <- plt2 + aes_string(shape = 'Estimates', linetype = 'Estimates', colour = 'Estimates')
          plt2 <- plt2 + geom_line(size = 0.8)
          plt2 <- plt2 + xlab(parse(text=x$xname[i]))
          plt2 <- plt2 + theme_bw()
          plt2 <- plt2 + theme(legend.position = "top")
          plt2 <- plt2 + ylab('Partial Residuals')
          plt2 <- plt2 + scale_linetype_manual(values = c("dotdash", "dotdash", "solid", "solid"))
        }
        grid.arrange(plt1, plt2, nrow = 2)
      }
    }
  } else {
    if (nfun == 1) {
      if (x$model == 'gbsar') {
        if (x$family == "bernoulli") {
          if (x$link == 'probit') {
            fxl = pnorm(fxobsl + rep(median(x$wbeta$upper), nobs))
            fxm = pnorm(fxobsm + rep(median(x$wbeta$mean), nobs))
            fxu = pnorm(fxobsu + rep(median(x$wbeta$lower), nobs))
          } else {
            logit <- function(xx) 1 / (1 + exp(-xx))
            fxl = logit(fxobsl + rep(median(x$wbeta$upper), nobs))
            fxm = logit(fxobsm + rep(median(x$wbeta$mean), nobs))
            fxu = logit(fxobsu + rep(median(x$wbeta$lower), nobs))
          }
          if (!ask) dev.new()
          plot(x = xgrid[, 1], y = fxgridm[, 1], main = '', pch = NA, ylim = range(x$fxgrid), 
               xlab = x$xname[1], ylab = paste('f(', x$xname, ')', sep = ''))
          polygon(x = c(xgrid[, 1], rev(xgrid[, 1])), 
                  y = c(fxgridl[, 1], rev(fxgridu[, 1])), col = 'gray70', lty = 2)
          lines(x = xgrid[, 1], y = fxgridm[, 1], lwd = 2, col = 2)
          if (!ask) dev.new()
          o = order(xobs[, 1])
          plot(x = xobs[, 1], y = yobs, pch = NA, xlab = x$xname[1], main = '',
               ylab = paste('P(', x$yname, ')', sep = ''), ylim = c(0, 1))
          polygon(x = c(xobs[o, 1], rev(xobs[o, 1])), 
                  y = c(fxl[o], rev(fxu[o])), col = 'gray70', lty = 2)
          lines(x = xobs[o, 1], y = fxm[o], lwd = 2, col = 2)
          rug(x = xobs[yobs == 0, 1], side = 1)
          rug(x = xobs[yobs == 1, 1], side = 3)
        } else {
          if (!ask) dev.new()
          plot(x = xgrid[, 1], y = fxgridm[, 1], main = '', pch = NA, ylim = range(x$fxgrid), 
               xlab = x$xname[1], ylab = paste('f(', x$xname, ')', sep = ''))
          polygon(x = c(xgrid[, 1], rev(xgrid[, 1])), 
                  y = c(fxgridl[, 1], rev(fxgridu[, 1])), col = 'gray70', lty = 2)
          lines(x = xgrid[, 1], y = fxgridm[, 1], lwd = 2, col = 2)
          if (!ask) dev.new()
          fxl = exp(fxobsl + rep(median(x$wbeta$upper), nobs))
          fxm = exp(fxobsm + rep(median(x$wbeta$mean), nobs))
          fxu = exp(fxobsu + rep(median(x$wbeta$lower), nobs))
          o = order(xobs[, 1])
          plot(x = xobs[, 1], y = yobs, pch = NA, xlab = x$xname[1], main = '',
               ylab = x$yname[1], ylim = range(c(yobs, fxl, fxu)))
          polygon(x = c(xobs[o, 1], rev(xobs[o, 1])), 
                  y = c(fxl[o], rev(fxu[o])), col = 'gray70', lty = 2)
          points(x = xobs[, 1], y = yobs, lwd = 2, pch = 1)
          lines(x = xobs[o, 1], y = fxm[o], lwd = 2, col = 2)
        }
      } else {
        resid <- yobs - wbm
        plot(xobs, resid, pch = NA, ylim = range(c(resid, fxgridu, fxgridl)),
             xlab = x$xname, ylab = 'Parametric Residuals', main = '')
        polygon(x = c(xgrid, rev(xgrid)),
                y = c(fxgridl, rev(fxgridu)), col = 'gray70', lty = 2)
        points(xobs, resid, lwd = 2)
        lines(xgrid, fxgridm, lwd = 3, lty = 1, col = 2)
      }
    } else {
      for (i in 1:nfun) {
        if (!ask) dev.new()
        plot(x = xgrid[, i], y = fxgridm[, i], pch = NA, 
             ylim = range(c(fxgridl[, i], fxgridu[, i])), main = '', 
             xlab = x$xname[i], ylab = paste('f(', x$xname[i], ')', sep = ''))
        polygon(x = c(xgrid[, i], rev(xgrid[, i])), 
                y = c(fxgridl[, i], rev(fxgridu[, i])), col = 'gray70', lty = 2)
        lines(x = xgrid[, i], y = fxgridm[, i], lwd = 2, col = 2)
        if (!ask) dev.new()
        if (x$model == 'gbsar') {
          if (x$family == 'bernoulli') {
            if (x$link == 'probit') {
              fxl = pnorm(fxobsl[, i])
              fxm = pnorm(fxobsm[, i])
              fxu = pnorm(fxobsu[, i])
            } else {
              logit <- function(xx) 1 / (1 + exp(-xx))
              fxl = logit(fxobsl[, i])
              fxm = logit(fxobsm[, i])
              fxu = logit(fxobsu[, i])
            }
            o = order(xobs[, i])
            plot(x = xobs[o, i], y = fxm[o], pch = NA, ylim = c(0, 1), main = '', 
                 xlab = x$xname[i], ylab = paste('P(', x$yname, ')', sep = ''))
            polygon(x = c(xobs[o, i], rev(xobs[o, i])), 
                    y = c(fxl[o], rev(fxu[o])), col = 'gray70', lty = 2)
            lines(x = xobs[o, i], y = fxm[o], lwd = 2, col = 2)
            rug(xobs[yobs == 0, i], side = 1)
            rug(xobs[yobs == 1, i], side = 3)
          } else {
            fxl = exp(fxobsl[, i])
            fxm = exp(fxobsm[, i])
            fxu = exp(fxobsu[, i])
            o = order(xobs[, i])
            plot(x = xobs[o, i], y = fxm[o], pch = NA, ylim = range(c(yobs, fxl, fxu)), main = '', 
                 xlab = x$xname[i], ylab = x$yname)
            polygon(x = c(xobs[o, i], rev(xobs[o, i])), 
                    y = c(fxl[o], rev(fxu[o])), col = 'gray70', lty = 2)
            points(x = xobs[, i], y = yobs, lwd = 2)
            lines(x = xobs[o, i], y = fxm[o], lwd = 2, col = 2)
          }
        } else {
          resid = yobs - wbm - rowSums(fxobsm[, -i, drop = FALSE])
          plot(x = xobs[, i], y = resid, pch = NA, ylim = range(c(resid, fxobsl[, i], fxobsu[, i])), 
               main = '', xlab = x$xname[i], ylab = 'Partial Residuals')
          polygon(x = c(xgrid[, i], rev(xgrid[, i])), 
                  y = c(fxgridl[, i], rev(fxgridu[, i])), col = 'gray70', lty = 2)
          points(x = xobs[, i], y = yobs - wbm - rowSums(fxobsm[, -i, drop = FALSE]), lwd = 2)
          lines(x = xgrid[, i], y = fxgridm[, i], lwd = 2, col = 2)
        }
      }
    }
  }
}
