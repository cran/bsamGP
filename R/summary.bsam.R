"summary.bsam" <- function(object, ...) {
  out <- list()
  out$model <- object$model
  out$fmodel <- object$fmodel
  out$fpm <- object$fpm
  out$nbasis <- object$nbasis
  out$nfun <- object$nfun
  out$nExtremes <- object$nExtremes
  out$xmin <- object$xmin
  out$xmax <- object$xmax
  out$nobs <- object$n
  if (object$model == "bsaq") {
    out$p <- object$p
  }
  out$ndimw <- object$ndimw
  out$nblow <- object$mcmc$nblow
  out$smcmc <- object$mcmc$smcmc
  out$nskip <- object$mcmc$nskip
  out$nmcmc <- out$nblow + out$nskip * out$smcmc
  out$imodmet <- object$imodmet
  out$pmet <- object$pmet
  out$rsquarey <- object$rsquarey
  out$spmadeq <- object$spmadeq
  if (object$model != "gbsar") {
    if (object$spmadeq) {
      out$lmarg.lm <- object$lmarg.lm
    }
  }
  if (object$model == "gbsar") {
    out$family <- object$family
    out$link <- object$link
    if (object$family == "negative.binomial" || object$family == "poisson.gamma") {
      out$kappam <- object$post.est$kappam
      out$kappas <- object$post.est$kappas
    }
  }
  out$marglik <- object$marglik
  out$lilm <- object$lmarg.gd
  out$lilnr <- object$lmarg.nr
  out$betam <- object$post.est$betam
  out$betas <- object$post.est$betas
  out$wnames <- object$wnames
  out$sigmam <- object$post.est$sigmam
  out$sigmas <- object$post.est$sigmas
  out$alpham <- object$post.est$alpham
  out$alphas <- object$post.est$alphas
  out$psim <- object$post.est$psim
  out$psis <- object$post.est$psis
  out$omegam <- object$post.est$omegam
  out$omegas <- object$post.est$omegas
  out$taum <- object$post.est$taum
  out$taus <- object$post.est$taus
  out$gammam <- object$post.est$gammam
  out$gammas <- object$post.est$gammas
  out$zetam <- object$post.est$zetam
  out$zetas <- object$post.est$zetas
  out$thetam <- object$post.est$thetam
  out$thetas <- object$post.est$thetas
  out$iflagpsi <- object$prior$iflagpsi
  out$iflagprior <- object$prior$iflagprior
  class(out) <- 'summary.bsam'
  out
}
