"print.blm" <- function(x, ...) {
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  if (x$model == "blq") {
    cat("Bayesian Linear Quantile Regression \n")
    cat("\n")
    cat("Model:", "\n")
    cat("Y = w'beta + ALD(p;0,sigma^2)", "\n")
  }
  if (x$model == "gblr") {
    if (x$family == "bernoulli" && x$link == "probit") {
      cat("Bayesian Probit Regression (Albert-Chib) \n")
      cat("\n")
      cat("Model:", "\n")
      cat("Y ~ Ber(p), E(Y) = p", "\n")
      cat("p = F(w'beta), F : normal cdf", "\n")
    } else if (x$family == "bernoulli" && x$link == "logit") {
      if (x$algorithm == "KS") {
        cat("Bayesian Logit Regression (Kolmogorov-Smirnov) \n")
      } else {
        cat("Bayesian Logit Regression (Adaptive Metropolis) \n")
      }
      cat("\n")
      cat("Model:", "\n")
      cat("Y ~ Ber(p), E(Y) = p", "\n")
      cat("logit(p) = w'beta", "\n")
    }
    if (x$family == "poisson") {
      cat("Bayesian Poisson Regression \n")
      cat("\n")
      cat("Model:", "\n")
      cat("Y ~ Poi(lam), E(Y) = lam", "\n")
      cat("log(lam) = w'beta", "\n")
    }
    if (x$family == "negative.binomial") {
      cat("Bayesian Negative-Binomial Regression \n")
      cat("\n")
      cat("Model:", "\n")
      cat("Y ~ NegBin(lam,kappa)", "\n")
      cat("log(lam) = w'beta", "\n")
    }
    if (x$family == "poisson.gamma") {
      cat("Bayesian Negative-Binomial Regression \n")
      cat("\n")
      cat("Model:", "\n")
      cat("Y ~ Poi(lambda), lambda ~ Ga(kappa,kappa/mu)", "\n")
      cat("log(mu) = w'beta", "\n")
    }
  }
  if (x$model == "blr") {
    cat("Bayesian Linear Regression\n")
    cat("\n")
    cat("Model:", "\n")
    cat("Y = w'beta + N(0,sigma^2)", "\n")
  }
  cat("\n")
}