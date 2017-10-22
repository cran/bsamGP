"print.summary.bsamdpm" <- function(x, ...) {
  cat("\n")
  if (x$model == "bsaqdpm") {
    cat("Quantile of interest                 = ", x$p, "\n")
  }
  cat("Number of Cosine basis functions     = ", x$nbasis, "\n")
  cat("Number of observations               = ", x$nobs, "\n")
  cat("Number of covariates (no intercept)  = ", x$ndimw, "\n")
  cat("\n")
  cat("MCMC transition draws                = ", x$nblow, "\n")
  cat("MCMC draws saved for estimation      = ", x$smcmc, "\n")
  cat("Save every nskip draws               = ", x$nskip, "\n")
  cat("MCMC draws total                     = ", x$nmcmc, "\n")
  cat("\n")
  cat("-----", "\n")
  cat("\n")
  for (i in 1:x$nfun) {
    if (x$fmodel[i] > 1) {
      cat("Function = ", i, "\n")
      cat("Proportion of Theta Accepted after burn-in  = ", x$pmet[i], "\n")
      cat("\n")
    }
  }
  cat("Log Pseudo Marginal Likelihood (LPML)", "\n")
  cat("LPML Mukhopadhyay & Gelfand                 = ", round(x$lpml, 4), "\n")
  cat("\n")
  cat("R-Square                                    = ", round(x$rsquarey, 4), "\n")
  cat("\n")
  cat("Number of Clusters", "\n")
  print(table(x$nclass)/x$smcmc)
  cat("\n")
  cat("Total mass", "\n")
  print(summary(x$tmass))
  cat("\n")
  
  if (x$location) {
    if (x$ndimw >= 1) {
      cat("beta", "\n")
      sout <- c("PostM", "PostStd", "PostM/STD")
      bout <- cbind(x$betam, x$betas, (x$betam/x$betas))
      colnames(bout) <- sout
      rownames(bout) <- x$wnames
      print(bout)
      cat("\n")
    }
  } else {
    cat("beta", "\n")
    sout <- c("PostM", "PostStd", "PostM/STD")
    bout <- cbind(x$betam, x$betas, (x$betam/x$betas))
    colnames(bout) <- sout
    rownames(bout) <- x$wnames
    print(bout)
    cat("\n")
  }
  
  for (i in 1:x$nfun) {
    cat("-----", "\n")
    cat("\n")
    cat("Function = ", i, sep = "", "\n")
    if (x$fmodel[i] == 1) {
      cat("Unrestricted model", "\n")
    } else if (x$fmodel[i] == 2 & x$fpm[i] == 1) {
      cat("Increasing function", "\n")
    } else if (x$fmodel[i] == 2 & x$fpm[i] == -1) {
      cat("Decreasing function", "\n")
    } else if (x$fmodel[i] == 3 & x$fpm[i] == 1) {
      cat("Increasing convex function", "\n")
    } else if (x$fmodel[i] == 3 & x$fpm[i] == -1) {
      cat("Decreasing concave function", "\n")
    } else if (x$fmodel[i] == 4 & x$fpm[i] == 1) {
      cat("Increasing concave function", "\n")
    } else if (x$fmodel[i] == 4 & x$fpm[i] == -1) {
      cat("Decreasing convex function", "\n")
    } else if (x$fmodel[i] == 5 & x$fpm[i] == 1) {
      cat("Increasing S (convex to concave) function", "\n")
    } else if (x$fmodel[i] == 5 & x$fpm[i] == -1) {
      cat("Decreasing S (concave to convex) function", "\n")
    } else if (x$fmodel[i] == 6 & x$fpm[i] == 1) {
      cat("Increasing Rotated S (concave to convex) function", "\n")
    } else if (x$fmodel[i] == 6 & x$fpm[i] == -1) {
      cat("Decreasing Rotated S (convex to concave) function", "\n")
    } else if (x$fmodel[i] == 7 & x$fpm[i] == 1) {
      cat("Inverted U: increasing to omega then decreasing", "\n")
    } else if (x$fmodel[i] == 7 & x$fpm[i] == -1) {
      cat("U Shaped: decreasing to omega then increasing", "\n")
    }
    cat("\n")
    
    if (x$fmodel[i] > 2 & x$fmodel[i] < 7) {
      cat("Linear term alpha in x for constrainted f", "\n")
      cat("Posterior mean   of alpha = ", x$alpham[i], "\n")
      cat("Posterior stddev of alpha = ", x$alphas[i], "\n")
      cat("\n")
    }
    
    if ((x$fmodel[i] == 5) || (x$fmodel[i] == 6) || (x$fmodel[i] == 7)) {
      if (x$iflagpsi == 1) {
        cat("psi is slople of squish function", "\n")
        cat("Posterior mean   of psi = ", x$psim[i], "\n")
        cat("Posterior stddev of psi = ", x$psis[i], "\n")
      } else {
        cat("psi is slople of squish function, is constant", "\n")
        cat("Fixed psi = ", x$psim, "\n")
      }
      cat("\n")
      cat("omega is inflection point of squish function", "\n")
      cat("Posterior mean  omega = ", x$omegam[i], "\n")
      cat("Posterior stdev omega = ", x$omegas[i], "\n")
      cat("\n")
    }
    cat("theta_k ~ N(0,tau2*exp(-gamma*k))", "\n")
    cat("\n")
    if (x$fmodel[i] == 1) {
      if (x$iflagprior == 1) {
        cat("Lasso Smoothing prior with tau2 ~ Exp(lambda)", "\n")
      } else {
        cat("T-Smoothing prior with tau2 ~ IG", "\n")
      }
    }
    cat("Tau", "\n")
    cat("PostM   PostS", "\n")
    cat(c(x$taum[i], x$taus[i]), "\n")
    cat("\n")
    cat("Gamma", "\n")
    gamout <- c(x$gammam[i], x$gammas[i])
    names(gamout) <- c("PostM", "PostS")
    print(gamout)
    cat("\n")
    cat("Zeta = ln(tau2) - mean(k)*gamma", "\n")
    zetaout <- c(x$zetam[i], x$zetas[i])
    names(zetaout) <- c("PostM", "PostS")
    print(zetaout)
    
    cat("\n")
    cat("Cosine Basis weights theta", "\n")
    if (x$fmodel[i] == 1) {
      tnames <- paste("T", 1:x$nbasis, sep = "")
      bout <- cbind(x$thetam[2:(x$nbasis + 1), i], x$thetas[2:(x$nbasis + 1), i], 
                    (x$thetam[2:(x$nbasis + 1), i]/x$thetas[2:(x$nbasis + 1), i]))
    } else {
      tnames <- paste("T", 0:x$nbasis, sep = "")
      bout <- cbind(x$thetam[, i], x$thetas[, i], (x$thetam[, i]/x$thetas[, i]))
    }
    sout <- c("PostMean", "PostSTD", "Ratio")
    colnames(bout) <- sout
    rownames(bout) <- tnames
    print(round(bout, 5))
  }
  cat("\n")
}
