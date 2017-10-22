"summary.blm" <- function(object, ...) {
  out <- list()
  out$model <- object$model
  if (object$model == "blq") {
    out$p <- object$p
  }
  if (object$model == "gblr") {
    out$family <- object$family
    out$link <- object$link
  }
  out$n <- object$n
  out$ndimw <- object$ndimw
  out$nparw <- object$nparw
  out$nblow <- object$mcmc$nblow
  out$smcmc <- object$mcmc$smcmc
  out$nskip <- object$mcmc$nskip
  out$nmcmc <- out$nblow + out$nskip * out$smcmc
  out$lmarg <- object$lmarg
  out$betam <- object$post.est$betam
  out$betas <- object$post.est$betas
  out$wnames <- object$wnames
  if (object$model != "gblr") {
    out$sigmam <- object$post.est$sigmam
    out$sigmas <- object$post.est$sigmas
    out$rsquarey <- object$rsquarey
  }
  if (object$model == "gblr") {
    if (object$family == "negative.binomial" || object$family == "poisson.gamma") {
      out$kappam <- object$post.est$kappam
      out$kappas <- object$post.est$kappas
    }
  }
  out$marglik <- object$marglik
  class(out) <- "summary.blm"
  out
}
