"print.summary.blm" <- function(x, ...) {
  if (x$model == "blq") {
    cat("\n")
    cat("Quantile of interest                 = ", x$p, "\n")
  }
  cat("\n")
  cat("Number of observations               = ", x$n, "\n")
  cat("Number of covariates (no intercept)  = ", x$ndimw, "\n")
  cat("\n")
  cat("MCMC transition draws                = ", x$nblow, "\n")
  cat("MCMC draws saved for estimation      = ", x$smcmc, "\n")
  cat("Save every nskip draws               = ", x$nskip, "\n")
  cat("MCMC draws total                     = ", x$nmcmc, "\n")
  cat("\n")
  if (x$model != "gblr") {
    cat("R-Square                             = ", round(x$rsquarey, 4), "\n")
    cat("\n")
  }
  if (x$marglik) {
    cat("Log Integrated Likelihood", "\n")
    cat("LIL Parametric                       = ", round(x$lmarg, 4), "\n")
    cat("\n")
  }
  
  cat("beta", "\n")
  sout <- c("PostM", "PostStd", "PostM/STD")
  bout <- cbind(x$betam, x$betas, (x$betam/x$betas))
  colnames(bout) <- sout
  rownames(bout) <- x$wnames
  print(bout)
  cat("\n")
  if (x$model == "gblr") {
    if (x$family == "negative.binomial" || x$family == "poisson.gamma") {
      cat("kappa", "\n")
      cat("PostM kappa\t= ", x$kappam, "\n")
      cat("PostS kappa\t= ", x$kappas, "\n")
      cat("\n")
    }
  }
  if (x$model != "gblr") {
    cat("sigma", "\n")
    cat("PostM sigma\t= ", x$sigmam, "\n")
    cat("PostS sigma\t= ", x$sigmas, "\n")
    cat("\n")
  }
}