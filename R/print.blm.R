"print.blm" <- function(x, ...) {
  nobs <- x$n
  if (x$model == "blq") 
    p <- x$p
  if (x$model == "gblr") {
    family <- x$family
    link <- x$link
  }
  ndimw <- x$ndimw
  nparw <- x$nparw
  nblow <- x$mcmc$nblow
  smcmc <- x$mcmc$smcmc
  nskip <- x$mcmc$nskip
  nmcmc <- nblow + nskip * smcmc
  lmarg.lm <- x$lmarg
  betam <- x$post.est$betam
  betas <- x$post.est$betas
  wnames <- x$wnames
  if (x$model != "gblr") {
    sigmam <- x$post.est$sigmam
    sigmas <- x$post.est$sigmas
    rsquarey <- x$rsquarey
  }
  if (x$model == "gblr") {
    if (family == "negative.binomial" || family == "poisson.gamma") {
      kappam <- x$post.est$kappam
      kappas <- x$post.est$kappas
    }
  }
  
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  if (x$model == "blq") {
    cat("Bayesian Linear Quantile Regression \n")
    cat("\n")
    cat("\tModel:", "\n")
    cat("  \tY = w'beta + ALD(p;0,sigma^2)", "\n")
  }
  if (x$model == "gblr") {
    if (family == "bernoulli" && link == "probit") {
      cat("Bayesian Probit Regression (Albert-Chib) \n")
      cat("\n")
      cat("Model:", "\n")
      cat("     Y ~ Ber(p), E(Y) = p", "\n")
      cat("     p = F(w'beta), F : normal cdf", "\n")
    } else if (family == "bernoulli" && link == "logit") {
      cat("Bayesian Logit Regression (Kolmogorov-Smirnov) \n")
      cat("\n")
      cat("Model:", "\n")
      cat("     Y ~ Ber(p), E(Y) = p", "\n")
      cat("     logit(p) = w'beta", "\n")
    }
  }
  if (x$model == "blr") {
    cat("Bayesian Linear Regression\n")
    cat("\n")
    cat("\tModel:", "\n")
    cat("  \tY = w'beta + N(0,sigma^2)", "\n")
  }
  cat("\t\t", "\n")
  cat("**********************************************************************************", "\n")
  cat("\n")
  if (x$model == "blq") {
    cat("Quantile of interest                 = ", p, "\n")
  }
  cat("Number of observations               = ", nobs, "\n")
  cat("Number of covariates (no intercept)  = ", ndimw, "\n")
  cat("\n")
  cat("MCMC transition draws                = ", nblow, "\n")
  cat("MCMC draws saved for estimation      = ", smcmc, "\n")
  cat("Save every nskip draws               = ", nskip, "\n")
  cat("MCMC draws total                     = ", nmcmc, "\n")
  cat("\n")
  if (x$model != "gblr") {
    cat("R-Square                             = ", round(rsquarey, 4), "\n")
    cat("\n")
  }
  if (x$marglik) {
    cat("Log Integrated Likelihood", "\n")
    cat("LIL Parametric                       = ", round(lmarg.lm, 4), "\n")
    cat("\n")
  }
  
  cat("beta", "\n")
  cat("CNST is Y intercept", "\n")
  sout <- c("PostM", "PostStd", "PostM/STD")
  bout <- cbind(betam, betas, (betam/betas))
  colnames(bout) <- sout
  rownames(bout) <- wnames
  print(bout)
  cat("\n")
  if (x$model == "gblr") {
    if (family == "negative.binomial" || family == "poisson.gamma") {
      cat("kappa", "\n")
      cat("PostM kappa\t= ", kappam, "\n")
      cat("PostS kappa\t= ", kappas, "\n")
      cat("\n")
    }
  }
  if (x$model != "gblr") {
    cat("sigma", "\n")
    cat("PostM sigma\t= ", sigmam, "\n")
    cat("PostS sigma\t= ", sigmas, "\n")
    cat("\n")
  }
  cat("**********************************************************************************", "\n")
  cat("\n")
}