'summary.bsad' <- function(object, ...)
{
  out <- list()
  out$kappaloop <- object$mcmc$kappaloop
  out$nblow <- object$mcmc$nblow
  out$nskip <- object$mcmc$nskip
  out$smcmc <- object$mcmc$smcmc
  out$nobs <- object$nobs
  out$xmin <- object$xmin
  out$xmax <- object$xmax
  out$nint <- object$nint
  out$MaxNCos <- object$MaxNCos
  out$parametric <- object$parametric
  out$PostProbs <- object$PostProbs
  out$marginal.likelihood <- object$marginal.likelihood
  out$lmarg <- object$lmarg
  out$betaParm <- object$post.est$betaParm
  out$betaPars <- object$post.est$betaPars
  out$betaParql <- apply(object$mcmc.draws$betaPar, 2, function(x) quantile(x,probs = 0.025))
  out$betaParq2 <- apply(object$mcmc.draws$betaPar, 2, function(x) quantile(x,probs = 0.5))
  out$betaParqu <- apply(object$mcmc.draws$betaPar, 2, function(x) quantile(x,probs = 0.975))
  out$betam <- object$post.est$betam
  out$betas <- object$post.est$betas
  out$betaql <- apply(object$mcmc.draws$beta,2,function(x) quantile(x,probs=0.025))
  out$betaq2 <- apply(object$mcmc.draws$beta,2,function(x) quantile(x,probs=0.5))
  out$betaqu <- apply(object$mcmc.draws$beta,2,function(x) quantile(x,probs=0.975))
  out$kappam <- object$post.est$kappam
  out$kappas <- object$post.est$kappas
  out$kappaq <- quantile(object$mcmc.draws$kappa,probs=c(0.025,0.5,0.975))
  out$taum <- object$post.est$taum
  out$taus <- object$post.est$taus
  out$tauq <- quantile(object$mcmc.draws$tau,probs=c(0.025,0.5,0.975))
  out$gammam <- object$post.est$gamm
  out$gammas <- object$post.est$gams
  out$gammaq <- quantile(object$mcmc.draws$gam,probs=c(0.025,0.5,0.975))
  out$betaMaxKappam <- object$post.est$betaMaxKappam
  out$betaMaxKappas <- object$post.est$betaMaxKappas
  out$betaMaxKappaql <- apply(object$mcmc.draws$betaMaxKappa,2,function(x) quantile(x,probs=0.025))
  out$betaMaxKappaq2 <- apply(object$mcmc.draws$betaMaxKappa,2,function(x) quantile(x,probs=0.5))
  out$betaMaxKappaqu <- apply(object$mcmc.draws$betaMaxKappa,2,function(x) quantile(x,probs=0.975))
  out$tauMaxKappam <- object$post.est$tauMaxKappam
  out$tauMaxKappas <- object$post.est$tauMaxKappas
  out$tauMaxKappaq <- quantile(object$mcmc.draws$tauMaxKappa, probs=c(0.025,0.5,0.975))
  out$gamMaxKappam <- object$post.est$gamMaxKappam
  out$gamMaxKappas <- object$post.est$gamMaxKappas
  out$gamMaxKappaq <- quantile(object$mcmc.draws$gamMaxKappa, probs=c(0.025,0.5,0.975))
  class(out) <- 'summary.bsad'
  out
}
