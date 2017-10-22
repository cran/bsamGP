"summary.bsamdpm" <- function(object, ...) {
  out <- list()
  out$model <- object$model
  out$fmodel <- object$fmodel
  out$fpm <- object$fpm
  out$nbasis <- object$nbasis
  out$ndimw <- object$ndimw
  out$nfun <- object$nfun
  out$xmin <- object$xmin
  out$xmax <- object$xmax
  out$nobs <- object$n
  if (object$model == "bsaqdpm") {
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
  out$betam <- object$post.est$betam
  out$betas <- object$post.est$betas
  out$wnames <- object$wnames
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
  out$nclass <- object$dpm.draws$nclass
  out$tmass <- object$dpm.draws$tmass
  out$lpml <- object$lpml
  out$location <- object$location
  class(out) <- 'summary.bsamdpm'
  out
}
