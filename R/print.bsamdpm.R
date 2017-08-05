"print.bsamdpm" <- function(x, ...) {
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
  if (x$model == "bsaqdpm")
    p <- x$p
  ndimw <- x$ndimw
  nblow <- x$mcmc$nblow
  smcmc <- x$mcmc$smcmc
  nskip <- x$mcmc$nskip
  nmcmc <- nblow + nskip * smcmc
  imodmet <- x$imodmet
  pmet <- x$pmet
  rsquarey <- x$rsquarey
  betam <- x$post.est$betam
  betas <- x$post.est$betas
  wnames <- x$wnames
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
  nclass <- x$dpm.draws$nclass
  tmass <- x$dpm.draws$tmass
  lpml <- x$lpml

  if (x$model == "bsaqdpm") {
    cat("Bayesian Spectral Analysis Quantile Regression with Dirichlet process mixture errors \n")
    cat("\n")
    cat("\tModel:", "\n")
    cat("  \tY = w'beta + sum_k f_k(x_k) + int ALD(p;0,sigma^2)dG(sigma^2)", "\n")
    cat("  \tG ~ DP(M,G0), G0=Ga(1/sigma^2; r0/2,s0/2)", "\n")
  } else {
    cat("Bayesian Spectral Analysis Regression with Dirichlet process mixture errors \n")
    cat("\n")
    cat("\tModel:", "\n")
    if (x$location) {
      cat("        Y = w'beta + sum_k f_k(x_k) + int N(mu,sigma^2)dG(mu,sigma^2)", "\n")
      cat("        G ~ DP(M,G0),                                 ", "\n")
      cat("\t\tG0=N(mu;0,kappa/sigma^2)Ga(1/sigma^2; r0/2,s0/2)", "\n")
    } else {
      cat("  \tY = w'beta + sum_k f_k(x_k) + int N(0,sigma^2)dG(sigma^2)", "\n")
      cat("  \tG ~ DP(M,G0), G0=Ga(1/sigma^2; r0/2,s0/2)", "\n")
    }
  }
  cat("\t\txmin < x < xmax", "\n")
  cat("\t\t  Normalize f: int_a^b f(x)dx = 0", "\n")
  cat("\t\tbeta has the Intercept", "\n")
  cat("\tDefine", "\n")
  cat("\t\tZ(x) \t= sum_j=0^J theta_j*phi_j(x)", "\n")
  cat("\t\tg^a(x) \t= int_0^x |Z(s)|^2 dx", "\n")
  cat("\t\t\t \t= theta' Phi^a(x)*theta", "\n")
  cat("\t\twhere Phi^a(x) is matrix with ", "\n")
  cat("\t\tphi^a_{j,k}(x) = int_0^x phi_j(s) phi_k(s) ds in (j,k) ", "\n")
  cat("\t\tg^b(x)\t= int_0^x g^a(s) ds", "\n")
  cat("\t\t\t\t= theta' Phi^b(x)*theta", "\n")
  cat("\t\twhere Phi^b(x) is a matrix with", "\n")
  cat("\t\tphi^b_{j,k}(x) = int_0^x phi^a_{j,k}(s) ds", "\n")
  cat("\t\t", "\n")
  cat("*************************************************************************************", "\n")
  cat("\n")
  if (x$model == "bsaqdpm") {
    cat("Number of BSAQ functions = ", nfun, "\n")
  } else {
    cat("Number of BSAR functions = ", nfun, "\n")
  }
  cat(" Model for f", "\n")
  for (i in 1:nfun) {
    if (fmodel[i] == 1) {
      cat("\t\tNo shape restriction", "\n")
      cat("\t\t\tf(x) = Z(x)  ", "\n")
      cat("\t\tFree f does not have theta_0", "\n")
    } else {
      cat("Removed theta_0^2 from f but kept theta_0*theta_j", "\n")
    }

    if (fmodel[i] == 2) {
      if (fpm[i] == 1) {
        cat("\t\tIncreasing f", "\n")
        cat("\t\t\tf(x) = g^a(x) ", "\n")
      } else {
        cat("\t\tDecreasing f", "\n")
        cat("\t\t\tf(x) = -g^a(x) ", "\n")
      }
    }

    if (fmodel[i] == 3) {
      if (fpm[i] == 1) {
        cat("\t\tIncreasing convex f", "\n")
        cat("\t\t\tf(x) = g^b(x) +  alpha*(x-xmid)", "\n")
        cat("\t\t\talpha > 0", "\n")
      } else {
        cat("\t\tDecreasing, concave f", "\n")
        cat("\t\t\tf(x) = -g^b(x) +  alpha*(x-xmid)", "\n")
        cat("\t\t\talpha < 0", "\n")
      }
    }

    if (fmodel[i] == 4) {
      if (fpm[i] == 1) {
        cat("\t\tIncreasing concave f", "\n")
        cat("\t\t\tf(x) = -g^b(a+b-x) +  alpha*(x-xmid) ", "\n")
        cat("\t\t\talpha > 0", "\n")
      } else {
        cat("\t\tDecreasing convex f", "\n")
        cat("\t\t\tf(x) = g^b(a+b-x) +  alpha*(x-xmid)", "\n")
        cat("\t\t\talpha < 0", "\n")
      }
    }

    if (fmodel[i] == 5) {
      if (fpm[i] == 1) {
        cat("\t\tIncreasing, S shaped", "\n")
        cat("\t\t\tf(x) = int_a^x int_a^s Z(t)^2 h(t)dt ds +  alpha*(x-xmid) - xi*(x-xmin)", "\n")
        cat("\t\t\th(t) = {1-exp[psi*(x-omega)]}/{1+exp[psi*(x-omega)]} ", "\n")
        cat("\t\t\tpsi  > 0, a < omega < b, alpha > 0", "\n")
        cat("\t\t\txi = min(0,f'(x)) to make sure that f' > 0", "\n")
      } else {
        cat("\t\tDecreasing, S shaped", "\n")
        cat("\t\t\tf(x) = -int_a^x int_a^s Z(t)^2 h(t)dt ds +  alpha*(x-xmid) - xi*(x-xmin)", "\n")
        cat("\t\t\th(t) = {1-exp[psi*(x-omega)]}/{1+exp[psi*(x-omega)]} ", "\n")
        cat("\t\t\txi   = max(0,f'(x)) to make sure that f' < 0", "\n")
        cat("\t\t\tpsi  > 0, a < omega < b, alpha < 0", "\n")
      }
    }

    if (fmodel[i] == 6) {
      if (fpm[i] == 1) {
        cat("\t\tIncreasing concave-to-convex", "\n")
        cat("\t\t\tConcave before omega and convex after omega", "\n")
        cat("\t\t\tf(x) = int_a^x int_a^s Z(t)^2 h2(t)dt ds +  alpha*(x-xmid) - xi*(x-xmin)", "\n")
        cat("\t\t\th2(t) = {exp[psi*(x-omega)]-1}/{exp[psi*(x-omega)]+1}", "\n")
        cat("\t\t\tpsi  > 0, a < omega < b, alpha > 0", "\n")
        cat("\t\t\txi = min(0,f'(x)) to make sure that f' > 0", "\n")
      } else {
        cat("\t\tDecreasing convex-to-concave", "\n")
        cat("\t\t\tConcave before omega and convex after omega", "\n")
        cat("\t\t\tf(x) = int_a^x int_a^s Z(t)^2 h2(t)dt ds +  alpha*(x-xmid) - xi*(x-xmin)", "\n")
        cat("\t\t\th2(t) = {exp[psi*(x-omega)]-1}/{exp[psi*(x-omega)]+1} ", "\n")
        cat("\t\t\tpsi  > 0, a < omega < b, alpha > 0", "\n")
        cat("\t\t\txi = min(0,f'(x)) to make sure that f' > 0", "\n")
      }
    }

    if (fmodel[i] == 7) {
      if (fpm[i] == 1) {
        cat("\t\tConcave increasing to decreasing or inverted U shaped", "\n")
        cat("\t\t\tIncreasing before omega and decreasing after omega", "\n")
        cat("\t\t\tf(x) \t= int_0^x Z^2(s)h(s) ds - zeta + alpha_0 ", "\n")
        cat("\t\t\tzeta \t= min(0,min int_0^x Z^2(s)h(s)ds)", "\n")
        cat("\t\t\th(x) \t= {1- exp[psi*(x-omega)]}/{exp[psi*(x-omega)]+1} ", "\n")
      } else {
        cat("\t\tConvex decreasing to increasing or U shaped", "\n")
        cat("\t\t\tf(x) \t= -int_0^x Z^2(s)h(s) ds + zeta + alpha_0 ", "\n")
        cat("\t\t\tzeta \t= min(0,min int_0^x Z^2(s)h(s)ds)", "\n")
        cat("\t\t\th(x) \t= {1 - exp[psi*(x-omega)]}/{exp[psi*(x-omega)]+1}", "\n")
      }
    }

    if (fmodel[i] > 2 & fmodel[i] < 7) {
      cat("Model for f includes linear term alpha", "\n")
    }
  }
  cat("\n")
  cat("*************************************************************************************", "\n")
  cat("", "\n")
  cat("\tCosine basis ", "\n")
  cat("\txrange   = xmax-xmin", "\n")
  cat("\tphi_0(x) = 1/sqrt(xmax-xmin) for xmin < x < xmax", "\n")
  cat("\tphi_k(x) = sqrt(2/(xmax-xmin))*cos(pi*k*(x-xmin)/(xmax-xmin))", "\n")
  cat("", "\n")
  cat("\tScale invariant priors:", "\n")
  cat("\tbeta \t\t~ N(b0,B0)", "\n")
  if (sum(fmodel > 2) > 0)
    cat("\talpha\t\t~ N(m0,v0)I(delta*alpha>0)", "\n")
  if (sum(fmodel == 1) > 0)
    cat("\ttheta_0 \t~ N(0,v0)", "\n")
  if (sum(fmodel >= 2) > 0)
    cat("\ttheta_0\t\t~ N(0,v0)I(theta_0 > 0)", "\n")
  if (sum(fmodel >= 1) > 0)
    cat("\ttheta_k \t~ N(0,tau^2*exp(-gamma*k)", "\n")
  cat("", "\n")
  cat("\tSmoothing parameters tau and gamma", "\n")
  cat("\tChoice of two priors for tau2", "\n")
  if (x$prior$iflagprior == 0)
    cat("\t\t1.) T Smoother: \ttau2 ~ IG(r0/2,s0/2)", "\n")
  if (x$prior$iflagprior == 1)
    cat("\t\t2.) Lasso Smoother:\ttau2 ~ Exp(u0)", "\n")
  cat("\tgamma ~ Exp(w0)", "\n")
  cat("\tNote: posterior of tau and gamma have banana contours", "\n")
  cat("\tzeta = ln(tau2) - kbar*gamma is less dependent with gamma or log(gamma)", "\n")
  cat("", "\n")
  if ((fmodel[i] == 5) || (fmodel[i] == 6) || (fmodel[i] == 7)) {
    cat("\tS models uses squish function (reparameterization of hyperbolic tangent)", "\n")
    cat("\tthat depens on slope psi and location omega (inflection point for S).", "\n")
    cat("\tpsi \t~ N(m0,v0)I(psi > 0)  Truncated Normal", "\n")
    cat("  omega \t~ N(m0,v0)I(xmin < omega < xmax) and m0 = (xmin+xmax)/2", "\n")
    cat("\t", "\n")
  }
  cat("*************************************************************************************", "\n")
  cat("\n")
  if (x$model == "bsaqdpm") {
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
  cat("-------------------------------------------------------------------------------------", "\n")
  cat("\n")
  for (i in 1:nfun) {
    if (fmodel[i] > 1) {
      cat("Function = ", i, "\n")
      cat("Proportion of Theta Accepted after burn-in  = ", pmet[i], "\n")
      cat("\n")
    }
  }
  cat("Log Pseudo Marginal Likelihood (LPML)", "\n")
  cat("LPML Mukhopadhyay & Gelfand                 = ", round(lpml, 4), "\n")
  cat("\n")
  cat("R-Square                                    = ", round(rsquarey, 4), "\n")
  cat("\n")
  cat("Number of Clusters", "\n")
  print(table(nclass)/smcmc)
  cat("\n")
  cat("Total mass", "\n")
  print(summary(tmass))
  cat("\n")

  if (x$location) {
    if (x$ndimw >= 1) {
      cat("beta", "\n")
      sout <- c("PostM", "PostStd", "PostM/STD")
      bout <- cbind(betam, betas, (betam/betas))
      colnames(bout) <- sout
      rownames(bout) <- wnames
      print(bout)
      cat("\n")
    }
  } else {
    cat("beta", "\n")
    cat("CNST is Y intercept", "\n")
    sout <- c("PostM", "PostStd", "PostM/STD")
    bout <- cbind(betam, betas, (betam/betas))
    colnames(bout) <- sout
    rownames(bout) <- wnames
    print(bout)
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
    cat("theta_k ~ N(0,tau2*exp(-gamma*k))", "\n")
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
  cat("\n")
}
