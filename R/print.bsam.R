"print.bsam" <- function(x, ...) {
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")

  fmodel <- x$fmodel
  fpm <- x$fpm
  nbasis <- x$nbasis
  nfun <- x$nfun
  xmin <- x$xmin
  xmax <- x$xmax
  nobs <- x$n
  if (x$model == "bsaq")
    p <- x$p
  ndimw <- x$ndimw
  nblow <- x$mcmc$nblow
  smcmc <- x$mcmc$smcmc
  nskip <- x$mcmc$nskip
  nmcmc <- nblow + nskip * smcmc
  imodmet <- x$imodmet
  pmet <- x$pmet
  rsquarey <- x$rsquarey
  if (x$model != "gbsar") {
    if (x$spmadeq)
      lmarg.lm <- x$lmarg.lm
  }
  if (x$model == "gbsar") {
    family <- x$family
    link <- x$link
    if (family == "negative.binomial" || family == "poisson.gamma") {
      kappam <- x$post.est$kappam
      kappas <- x$post.est$kappas
    }
  }
  lilm <- x$lmarg.gd
  lilnr <- x$lmarg.nr
  betam <- x$post.est$betam
  betas <- x$post.est$betas
  wnames <- x$wnames
  sigmam <- x$post.est$sigmam
  sigmas <- x$post.est$sigmas
  alpham <- x$post.est$alpham
  alphas <- x$post.est$alphas
  psim <- x$post.est$psim
  psis <- x$post.est$psis
  omegam <- x$post.est$omegam
  omegas <- x$post.est$omegas
  taum <- x$post.est$taum
  taus <- x$post.est$taus
  gammam <- x$post.est$gammam
  gammas <- x$post.est$gammas
  zetam <- x$post.est$zetam
  zetas <- x$post.est$zetas
  thetam <- x$post.est$thetam
  thetas <- x$post.est$thetas
  iflagpsi <- x$prior$iflagpsi
  iflagprior <- x$prior$iflagprior

  if (x$model == "bsaq") {
    cat("Bayesian Spectral Analysis Quantile Regression (BSAQ) \n")
    cat("\n")
    cat("Model:", "\n")
    cat("Y = w'beta + sum_k f_k(x_k) + ALD(p;0,sigma^2)", "\n")
  } else if (x$model == "gbsar") {
    if (family == "bernoulli" && link == "probit") {
      cat("Bayesian Spectral Analysis Probit Regression (Alber-Chib) \n")
      cat("\n")
      cat("Model:", "\n")
      cat("Y ~ Ber(p), E(Y) = p", "\n")
      cat("p = F(w'beta + sum_k f_k(x_k)), F : normal cdf", "\n")
    }
    if (family == "bernoulli" && link == "logit") {
      cat("Bayesian Spectral Analysis Logistic Regression (Kolmogorov-Smirnov) \n")
      cat("\n")
      cat("Model:", "\n")
      cat("Y ~ Ber(p), E(Y) = p", "\n")
      cat("logit(p) = w'beta + sum_k f_k(x_k)", "\n")
    }
    if (family == "poisson") {
      cat("Bayesian Spectral Analysis Poisson Regression \n")
      cat("\n")
      cat("Model:", "\n")
      cat("Y ~ Poi(lam), E(Y) = lam", "\n")
      cat("log(lam) = w'beta + sum_k f_k(x_k)", "\n")
    }
    if (family == "negative.binomial") {
      cat("Bayesian Spectral Analysis Negative-Binomial Regression \n")
      cat("\n")
      cat("Model:", "\n")
      cat("Y ~ NegBin(lam,kappa)", "\n")
      cat("log(lam) = w'beta + sum_k f_k(x_k)", "\n")
    }
    if (family == "poisson.gamma") {
      cat("Bayesian Spectral Analysis Negative-Binomial Regression \n")
      cat("\n")
      cat("Model:", "\n")
      cat("Y ~ Poi(lambda), lambda ~ Ga(kappa,kappa/mu)", "\n")
      cat("log(mu) = w'beta + sum_k f_k(x_k)", "\n")
    }
  } else {
    cat("Bayesian Spectral Analysis Regression (BSAR) \n")
    cat("\n")
    cat("Model:", "\n")
    cat("Y = w'beta + sum_k f_k(x_k) + N(0,sigma^2)", "\n")
  }
  cat("xmin < x < xmax", "\n")
  cat("Normalize f: int_a^b f(x)dx = 0", "\n")
  cat("beta has the Intercept", "\n")
  cat("\n")
  cat("Define", "\n")
  cat("Z(x) = sum_j=0^J theta_j*phi_j(x)", "\n")
  cat("g^a(x) = int_0^x |Z(s)|^2 dx", "\n")
  cat("       = theta' Phi^a(x)*theta", "\n")
  cat("where Phi^a(x) is matrix with ", "\n")
  cat("phi^a_{j,k}(x) = int_0^x phi_j(s) phi_k(s) ds in (j,k) ", "\n")
  cat("g^b(x)\t= int_0^x g^a(s) ds", "\n")
  cat("        = theta' Phi^b(x)*theta", "\n")
  cat("where Phi^b(x) is a matrix with", "\n")
  cat("phi^b_{j,k}(x) = int_0^x phi^a_{j,k}(s) ds", "\n")
  cat("\t\t", "\n")
  cat("*************************************************************************************", "\n")
  cat("\n")
  cat("Number of nonparametric components = ", nfun, "\n")
  cat("\n")
  cat("Model for f", "\n")
  for (i in 1:nfun) {
    if (fmodel[i] == 1) {
      cat("No shape restriction", "\n")
      cat("f(x) = Z(x)  ", "\n")
      cat("Free f does not have theta_0", "\n")
    } else {
      cat("Removed theta_0^2 from f but kept theta_0*theta_j", "\n")
    }

    if (fmodel[i] == 2) {
      if (fpm[i] == 1) {
        cat("Increasing f", "\n")
        cat("f(x) = g^a(x) ", "\n")
      } else {
        cat("Decreasing f", "\n")
        cat("f(x) = -g^a(x) ", "\n")
      }
    }

    if (fmodel[i] == 3) {
      if (fpm[i] == 1) {
        cat("Increasing convex f", "\n")
        cat("f(x) = g^b(x) +  alpha*(x-xmid)", "\n")
        cat("alpha > 0", "\n")
      } else {
        cat("Decreasing, concave f", "\n")
        cat("f(x) = -g^b(x) +  alpha*(x-xmid)", "\n")
        cat("alpha < 0", "\n")
      }
    }

    if (fmodel[i] == 4) {
      if (fpm[i] == 1) {
        cat("Increasing concave f", "\n")
        cat("f(x) = -g^b(a+b-x) +  alpha*(x-xmid) ", "\n")
        cat("alpha > 0", "\n")
      } else {
        cat("Decreasing convex f", "\n")
        cat("f(x) = g^b(a+b-x) +  alpha*(x-xmid)", "\n")
        cat("alpha < 0", "\n")
      }
    }

    if (fmodel[i] == 5) {
      if (fpm[i] == 1) {
        cat("Increasing, S shaped", "\n")
        cat("f(x) = int_a^x int_a^s Z(t)^2 h(t)dt ds +  alpha*(x-xmid) - xi*(x-xmin)", "\n")
        cat("h(t) = {1-exp[psi*(x-omega)]}/{1+exp[psi*(x-omega)]} ", "\n")
        cat("psi  > 0, a < omega < b, alpha > 0", "\n")
        cat("xi = min(0,f'(x)) to make sure that f' > 0", "\n")
      } else {
        cat("Decreasing, S shaped", "\n")
        cat("f(x) = -int_a^x int_a^s Z(t)^2 h(t)dt ds +  alpha*(x-xmid) - xi*(x-xmin)", "\n")
        cat("h(t) = {1-exp[psi*(x-omega)]}/{1+exp[psi*(x-omega)]} ", "\n")
        cat("xi   = max(0,f'(x)) to make sure that f' < 0", "\n")
        cat("psi  > 0, a < omega < b, alpha < 0", "\n")
      }
    }

    if (fmodel[i] == 6) {
      if (fpm[i] == 1) {
        cat("Increasing concave-to-convex", "\n")
        cat("Concave before omega and convex after omega", "\n")
        cat("f(x) = int_a^x int_a^s Z(t)^2 h2(t)dt ds +  alpha*(x-xmid) - xi*(x-xmin)", "\n")
        cat("h2(t) = {exp[psi*(x-omega)]-1}/{exp[psi*(x-omega)]+1}", "\n")
        cat("psi  > 0, a < omega < b, alpha > 0", "\n")
        cat("xi = min(0,f'(x)) to make sure that f' > 0", "\n")
      } else {
        cat("Decreasing convex-to-concave", "\n")
        cat("Concave before omega and convex after omega", "\n")
        cat("f(x) = int_a^x int_a^s Z(t)^2 h2(t)dt ds +  alpha*(x-xmid) - xi*(x-xmin)", "\n")
        cat("h2(t) = {exp[psi*(x-omega)]-1}/{exp[psi*(x-omega)]+1} ", "\n")
        cat("psi  > 0, a < omega < b, alpha > 0", "\n")
        cat("xi = min(0,f'(x)) to make sure that f' > 0", "\n")
      }
    }

    if (fmodel[i] == 7) {
      if (fpm[i] == 1) {
        cat("Concave increasing to decreasing or inverted U shaped", "\n")
        cat("Increasing before omega and decreasing after omega", "\n")
        cat("f(x) \t= int_0^x Z^2(s)h(s) ds - zeta + alpha_0 ", "\n")
        cat("zeta \t= min(0,min int_0^x Z^2(s)h(s)ds)", "\n")
        cat("h(x) \t= {1- exp[psi*(x-omega)]}/{exp[psi*(x-omega)]+1} ", "\n")
      } else {
        cat("Convex decreasing to increasing or U shaped", "\n")
        cat("f(x) \t= -int_0^x Z^2(s)h(s) ds + zeta + alpha_0 ", "\n")
        cat("zeta \t= min(0,min int_0^x Z^2(s)h(s)ds)", "\n")
        cat("h(x) \t= {1 - exp[psi*(x-omega)]}/{exp[psi*(x-omega)]+1}", "\n")
      }
    }

    if (fmodel[i] > 2 & fmodel[i] < 7) {
      cat("Model for f includes linear term alpha", "\n")
    }
    cat("\n")
    cat("*************************************************************************************", "\n")
  }
  cat("", "\n")
  cat("Cosine basis ", "\n")
  cat("xrange   = xmax-xmin", "\n")
  cat("phi_0(x) = 1/sqrt(xmax-xmin) for xmin < x < xmax", "\n")
  cat("phi_k(x) = sqrt(2/(xmax-xmin))*cos(pi*k*(x-xmin)/(xmax-xmin))", "\n")
  cat("", "\n")
  if (x$model == "gbsar") {
    cat("Priors:", "\n")
    cat("beta    ~ N(b0,B0)", "\n")
    if (sum(fmodel > 2) > 0)
      cat("alpha   ~ N(m0,v0)I(delta*alpha>0)", "\n")
    if (sum(fmodel == 1) > 0)
      cat("theta_0 ~ N(0,v0)", "\n")
    if (sum(fmodel >= 2) > 0)
      cat("theta_0 ~ N(0,v0)I(theta_0 > 0)", "\n")
    if (sum(fmodel == 1) > 0)
      cat("theta_k ~ N(0,tau^2*exp(-gamma*k)", "\n")
    if (sum(fmodel >= 2) > 0)
      cat("theta_k ~ N(0,tau^2*exp(-gamma*k)", "\n")
  } else {
    cat("Scale invariant priors:", "\n")
    cat("beta|sigma \t\t~ N(b0,sigma2*B0)", "\n")
    if (sum(fmodel > 2) > 0)
      cat("alpha|sigma\t\t~ N(m0,sigma2*v0)I(delta*alpha>0)", "\n")
    if (sum(fmodel == 1) > 0)
      cat("theta_0|sigma \t~ N(0,sigma2*v0)", "\n")
    if (sum(fmodel >= 2) > 0)
      cat("theta_0|sigma\t~ N(0,sigma*v0)I(theta_0 > 0)", "\n")
    if (sum(fmodel == 1) > 0)
      cat("theta_k|sigma \t~ N(0,sigma2*tau^2*exp(-gamma*k)", "\n")
    if (sum(fmodel >= 2) > 0)
      cat("theta_k|sigma \t~ N(0,sigma*tau^2*exp(-gamma*k)", "\n")
  }
  cat("", "\n")
  cat("Smoothing parameters tau and gamma", "\n")
  cat("Choice of two priors for tau2", "\n")
  if (x$prior$iflagprior == 0)
    cat("1.) T Smoother: tau2 ~ IG(r0/2,s0/2)", "\n")
  if (x$prior$iflagprior == 1)
    cat("2.) Lasso Smoother: tau2 ~ Exp(u0)", "\n")
  cat("gamma ~ Exp(w0)", "\n")
  cat("\n")
  cat("Note: posterior of tau and gamma have banana contours", "\n")
  cat("zeta = ln(tau2) - kbar*gamma is less dependent with gamma or log(gamma)", "\n")
  cat("", "\n")
  if ((fmodel[i] == 5) || (fmodel[i] == 6) || (fmodel[i] == 7)) {
    cat("S models uses squish function (reparameterization of hyperbolic tangent)", "\n")
    cat("that depens on slope psi and location omega (inflection point for S).", "\n")
    cat("psi \t~ N(m0,v0)I(psi > 0)  Truncated Normal", "\n")
    cat("omega \t~ N(m0,v0)I(xmin < omega < xmax) and m0 = (xmin+xmax)/2", "\n")
    cat("\t", "\n")
  }
  cat("*************************************************************************************", "\n")
  cat("\n")
  if (x$model == "bsaq") {
    cat("Quantile of interest                 = ", p, "\n")
  }
  cat("Number of Cosine basis functions     = ", nbasis, "\n")
  cat("Number of observations               = ", nobs, "\n")
  cat("Number of covariates (no intercept)  = ", ndimw, "\n")
  cat("\n")
  cat("MCMC transition draws                = ", nblow, "\n")
  cat("MCMC draws saved for estimation      = ", smcmc, "\n")
  cat("Save every nskip draws               = ", nskip, "\n")
  cat("MCMC draws total                     = ", nmcmc, "\n")
  cat("\n")
  for (i in 1:nfun) {
    if (fmodel[i] > 1) {
      cat("Function = ", i, "\n")
      cat("Proportion of Theta Accepted after burn-in    = ", pmet[i], "\n")
      cat("\n")
    }
  }
  if (x$model != "gbsar") {
    cat("R-Square                       = ", round(rsquarey, 4), "\n")
    cat("\n")
  }
  if (x$marglik) {
    cat("Log Integrated Likelihood", "\n")
    cat("LIL Gelfand & Dey              = ", round(lilm, 4), "\n")
    cat("LIL Newton & Raftery (biased)  = ", round(lilnr, 4), "\n")
    if (x$model != "gbsar") {
      if (x$spmadeq) {
        cat("LIL Parametric                 = ", round(lmarg.lm, 4), "\n")
        cat("\n")
        cat("H0: Parametric versus H1: Semiparametric", "\n")
        cat("Log Bayes Factor (BF[01])      = ", round(lilm - lmarg.lm, 4), "\n")
      }
    }
    cat("\n")
  }

  cat("beta", "\n")
  cat("const is Y intercept", "\n")
  sout <- c("PostM", "PostStd", "PostM/STD")
  bout <- cbind(betam, betas, (betam/betas))
  colnames(bout) <- sout
  rownames(bout) <- wnames
  print(bout)
  cat("\n")

  if (x$model == "gbsar") {
    if (family == "negative.binomial" || family == "poisson.gamma") {
      cat("kappa", "\n")
      cat("PostM kappa\t= ", kappam, "\n")
      cat("PostS kappa\t= ", kappas, "\n")
      cat("\n")
    }
  } else {
    cat("sigma", "\n")
    cat("PostM sigma\t= ", sigmam, "\n")
    cat("PostS sigma\t= ", sigmas, "\n")
    cat("\n")
  }

  for (i in 1:nfun) {
    cat("*************************************************************************************", "\n")
    cat("\n")
    cat("Function = ", i, sep = "", "\n")
    if (fmodel[i] == 1) {
      cat("Unrestricted model", "\n")
    } else if (fmodel[i] == 2 & fpm[i] == 1) {
      cat("Increasing function", "\n")
    } else if (fmodel[i] == 2 & fpm[i] == -1) {
      cat("Decreasing function", "\n")
    } else if (fmodel[i] == 3 & fpm[i] == 1) {
      cat("Increasing convex function", "\n")
    } else if (fmodel[i] == 3 & fpm[i] == -1) {
      cat("Decreasing concave function", "\n")
    } else if (fmodel[i] == 4 & fpm[i] == 1) {
      cat("Increasing concave function", "\n")
    } else if (fmodel[i] == 4 & fpm[i] == -1) {
      cat("Decreasing convex function", "\n")
    } else if (fmodel[i] == 5 & fpm[i] == 1) {
      cat("Increasing S (convex to concave) function", "\n")
    } else if (fmodel[i] == 5 & fpm[i] == -1) {
      cat("Decreasing S (concave to convex) function", "\n")
    } else if (fmodel[i] == 6 & fpm[i] == 1) {
      cat("Increasing Rotated S (concave to convex) function", "\n")
    } else if (fmodel[i] == 6 & fpm[i] == -1) {
      cat("Decreasing Rotated S (convex to concave) function", "\n")
    } else if (fmodel[i] == 7 & fpm[i] == 1) {
      cat("Inverted U: increasing to omega then decreasing", "\n")
    } else if (fmodel[i] == 7 & fpm[i] == -1) {
      cat("U Shaped: decreasing to omega then increasing", "\n")
    }
    cat("\n")

    if (fmodel[i] > 2 & fmodel[i] < 7) {
      cat("Linear term alpha in x for constrainted f", "\n")
      cat("Posterior mean   of alpha = ", alpham[i], "\n")
      cat("Posterior stddev of alpha = ", alphas[i], "\n")
      cat("\n")
    }

    if ((fmodel[i] == 5) || (fmodel[i] == 6) || (fmodel[i] == 7)) {
      if (iflagpsi == 1) {
        cat("psi is slople of squish function", "\n")
        cat("Posterior mean   of psi = ", psim[i], "\n")
        cat("Posterior stddev of psi = ", psis[i], "\n")
      } else {
        cat("psi is slople of squish function, is constant", "\n")
        cat("Fixed psi = ", psim, "\n")
      }
      cat("\n")
      cat("omega is inflection point of squish function", "\n")
      cat("Posterior mean  omega = ", omegam[i], "\n")
      cat("Posterior stdev omega = ", omegas[i], "\n")
      cat("\n")
    }
    if (x$model == "gbsar") {
      cat("theta_k ~ N(0,tau2*exp(-gamma*k))", "\n")
    } else {
      if (fmodel[i] == 1) {
        cat("theta_k ~ N(0,sigma2*tau2*exp(-gamma*k))", "\n")
      } else {
        cat("theta_k ~ N(0,sigma*tau2*exp(-gamma*k))", "\n")
      }
    }
    cat("\n")
    if (fmodel[i] == 1) {
      if (iflagprior == 1) {
        cat("Lasso Smoothing prior with tau2 ~ Exp(lambda)", "\n")
      } else {
        cat("T-Smoothing prior with tau2 ~ IG", "\n")
      }
    }
    cat("Tau", "\n")
    cat("PostM   PostS", "\n")
    cat(c(taum[i], taus[i]), "\n")
    cat("\n")
    cat("Gamma", "\n")
    gamout <- c(gammam[i], gammas[i])
    names(gamout) <- c("PostM", "PostS")
    print(gamout)
    cat("\n")
    cat("Zeta = ln(tau2) - mean(k)*gamma", "\n")
    zetaout <- c(zetam[i], zetas[i])
    names(zetaout) <- c("PostM", "PostS")
    print(zetaout)

    cat("\n")
    cat("Cosine Basis weights theta", "\n")
    if (fmodel[i] == 1) {
      tnames <- paste("T", 1:nbasis, sep = "")
      bout <- cbind(thetam[2:(nbasis + 1), i], thetas[2:(nbasis + 1), i], (thetam[2:(nbasis + 1), i]/thetas[2:(nbasis +
                                                                                                                 1), i]))
    } else {
      tnames <- paste("T", 0:nbasis, sep = "")
      bout <- cbind(thetam[, i], thetas[, i], (thetam[, i]/thetas[, i]))
    }
    sout <- c("PostMean", "PostSTD", "Ratio")
    colnames(bout) <- sout
    rownames(bout) <- tnames
    print(round(bout, 5))
  }
  cat("*************************************************************************************", "\n")
}
