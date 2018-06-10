"plot.fitted.bsad" <- function(x, ggplot2 = TRUE, legend.position = "top", nbins = 30, ...) {
  prob <- (1 - x$alpha) * 100
  HPD <- x$HPD
  par(...)
  if (x$parametric == "none") {
    if (ggplot2) {
      if (HPD) {
        datl <- data.frame(x = rep(x$xgrid, 3), dens = c(x$fsemi$upper, x$fsemi$mean, x$fsemi$lower),
                           Estimates = c(rep(paste(prob, "% HPD UCI (Semi)", sep = ""), x$nint),
                                         rep("Posterior Mean (Semi)", x$nint),
                                         rep(paste(prob, "% HPD LCI (Semi)", sep = ""), x$nint)))
      } else {
        datl <- data.frame(x = rep(x$xgrid, 3), dens = c(x$fsemi$upper, x$fsemi$mean, x$fsemi$lower),
                           Estimates = c(rep(paste(prob, "% Equal-tail UCI (Semi)", sep = ""), x$nint),
                                         rep("Posterior Mean (Semi)", x$nint),
                                         rep(paste(prob, "% Equal-tail LCI (Semi)", sep = ""), x$nint)))
      }
      datx = data.frame(x$x)

      plt1 <- ggplot(datl)
      plt1 <- plt1 + geom_histogram( data = datx, aes_string(x='x.x', y='..density..'), alpha=0.7,
                                     bins = nbins, col="black", fill= I("grey"), inherit.aes=FALSE)

      plt1 <- plt1 + aes_string(x = 'x', y = 'dens')
      plt1 <- plt1 + aes_string(group = 'Estimates')
      plt1 <- plt1 + aes_string(shape = 'Estimates', linetype = 'Estimates', colour = 'Estimates')
      plt1 <- plt1 + geom_line(size = 0.8)
      plt1 <- plt1 + xlab("")
      plt1 <- plt1 + ylab("Density")
      plt1 <- plt1 + theme_bw()
      plt1 <- plt1 + theme(legend.position = legend.position)
      plt1 <- plt1 + scale_linetype_manual(values = c("dotted", "dotted", "solid"))
      plt1
    } else {
      ymax <- max(c(x$fsemi$upper, hist(x$x, plot = F)$density))
      hist(x$x, prob = T, xlab = "", main = "", ylim = c(0, ymax), xlim = c(x$xmin, x$xmax), col = "gray95", nclass=nbins, ...)
      lines(x$xgrid, x$fsemi$mean, lwd = 2, lty = 1, col = "dodgerblue")
      lines(x$xgrid, x$fsemi$lower, lwd = 2, lty = 3, col = "tomato")
      lines(x$xgrid, x$fsemi$upper, lwd = 2, lty = 3, col = "seagreen3")
      rug(x$x, lwd = 2)
    }
  } else {
    if (ggplot2) {
      if (HPD) {
        datl <- data.frame(x = rep(x$xgrid, 4),
                           dens = c(x$fpar$mean, x$fsemi$upper, x$fsemi$mean, x$fsemi$lower),
                           Estimates = c(rep("Parametric", x$nint),
                                         rep(paste(prob, "% HPD UCI (Semi)", sep = ""), x$nint),
                                         rep("Posterior Mean (Semi)", x$nint),
                                         rep(paste(prob, "% HPD LCI (Semi)", sep = ""), x$nint)))
      } else {
        datl <- data.frame(x = rep(x$xgrid, 4),
                           dens = c(x$fpar$mean, x$fsemi$upper, x$fsemi$mean, x$fsemi$lower),
                           Estimates = c(rep("Parametric", x$nint),
                                         rep(paste(prob, "% Equal-tail UCI (Semi)", sep = ""), x$nint),
                                         rep("Posterior Mean (Semi)", x$nint),
                                         rep(paste(prob, "% Equal-tail LCI (Semi)", sep = ""), x$nint)))
      }
      datx = data.frame(x$x)
      plt1 <- ggplot(datl)
      plt1 <- plt1 + geom_histogram( data = datx, aes_string(x='x.x', y='..density..'), alpha=0.7,
                                     bins = nbins, col="black", fill= I("grey"), inherit.aes=FALSE)
      plt1 <- plt1 + aes_string(x = 'x', y = 'dens')
      plt1 <- plt1 + aes_string(group = 'Estimates')
      plt1 <- plt1 + aes_string(shape = 'Estimates', linetype = 'Estimates', colour = 'Estimates')
      plt1 <- plt1 + geom_line(size = 0.8)

      plt1 <- plt1 + xlab("")
      plt1 <- plt1 + ylab("Density")
      plt1 <- plt1 + theme_bw()
      plt1 <- plt1 + theme(legend.position = legend.position)
      plt1 <- plt1 + scale_linetype_manual(values = c("dotted", "dotted", "dashed", "solid"))
      plt1
    } else {
      ymax <- max(c(x$fsemi$upper, hist(x$x, plot = F)$density))
      hist(x$x, prob = T, xlab = "", main = "", col = "gray95", ylim = c(0, ymax), xlim = c(x$xmin, x$xmax), nclass = nbins, ...)
      lines(x$xgrid, x$fsemi$mean, lwd = 2, lty = 1, col = "magenta")
      lines(x$xgrid, x$fsemi$upper, lwd = 2, lty = 3, col = "seagreen3")
      lines(x$xgrid, x$fsemi$lower, lwd = 2, lty = 3, col = "tomato")
      lines(x$xgrid, x$fpar$mean, lwd = 2, lty = 4, col = "dodgerblue")
      rug(x$x, lwd = 2)
    }
  }
}
