"print.bsam" <- function(x, ...) {
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  if (x$model == "bsaq") {
    cat("Bayesian Spectral Analysis Quantile Regression (BSAQ) \n")
    cat("\n")
    cat("Model:", "\n")
    cat("Y = w'beta + sum_k f_k(x_k) + ALD(p;0,sigma^2)", "\n")
  } else if (x$model == "gbsar") {
    if (x$family == "bernoulli" && x$link == "probit") {
      cat("Bayesian Spectral Analysis Probit Regression (Alber-Chib) \n")
      cat("\n")
      cat("Model:", "\n")
      cat("Y ~ Ber(p), E(Y) = p", "\n")
      cat("p = F(w'beta + sum_k f_k(x_k)), F : normal cdf", "\n")
    }
    if (x$family == "bernoulli" && x$link == "logit") {
      if (x$algorithm == "KS") {
        cat("Bayesian Spectral Analysis Logistic Regression (Kolmogorov-Smirnov) \n")
      } else {
        cat("Bayesian Spectral Analysis Logistic Regression (Adaptive Metropolis) \n")
      }
      cat("\n")
      cat("Model:", "\n")
      cat("Y ~ Ber(p), E(Y) = p", "\n")
      cat("logit(p) = w'beta + sum_k f_k(x_k)", "\n")
    }
    if (x$family == "poisson") {
      cat("Bayesian Spectral Analysis Poisson Regression \n")
      cat("\n")
      cat("Model:", "\n")
      cat("Y ~ Poi(lam), E(Y) = lam", "\n")
      cat("log(lam) = w'beta + sum_k f_k(x_k)", "\n")
    }
    if (x$family == "negative.binomial") {
      cat("Bayesian Spectral Analysis Negative-Binomial Regression \n")
      cat("\n")
      cat("Model:", "\n")
      cat("Y ~ NegBin(lam,kappa)", "\n")
      cat("log(lam) = w'beta + sum_k f_k(x_k)", "\n")
    }
    if (x$family == "poisson.gamma") {
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
  cat("\n")
  cat("-----", "\n")
  cat("Number of nonparametric components = ", x$nfun, "\n")
  cat("\n")
  cat("Model for f", "\n")
  for (i in 1:x$nfun) {
    if (x$fmodel[i] == 1) {
      cat("No shape restriction", "\n")
      cat("f(x) = Z(x)  ", "\n")
      cat("Free f does not have theta_0", "\n")
    } else {
      cat("Removed theta_0^2 from f but kept theta_0*theta_j", "\n")
    }

    if (x$fmodel[i] == 2) {
      if (x$fpm[i] == 1) {
        cat("Increasing f", "\n")
        cat("f(x) = g^a(x) ", "\n")
      } else {
        cat("Decreasing f", "\n")
        cat("f(x) = -g^a(x) ", "\n")
      }
    }

    if (x$fmodel[i] == 3) {
      if (x$fpm[i] == 1) {
        cat("Increasing convex f", "\n")
        cat("f(x) = g^b(x) +  alpha*(x-xmid)", "\n")
        cat("alpha > 0", "\n")
      } else {
        cat("Decreasing, concave f", "\n")
        cat("f(x) = -g^b(x) +  alpha*(x-xmid)", "\n")
        cat("alpha < 0", "\n")
      }
    }

    if (x$fmodel[i] == 4) {
      if (x$fpm[i] == 1) {
        cat("Increasing concave f", "\n")
        cat("f(x) = -g^b(a+b-x) +  alpha*(x-xmid) ", "\n")
        cat("alpha > 0", "\n")
      } else {
        cat("Decreasing convex f", "\n")
        cat("f(x) = g^b(a+b-x) +  alpha*(x-xmid)", "\n")
        cat("alpha < 0", "\n")
      }
    }

    if (x$fmodel[i] == 5) {
      if (x$fpm[i] == 1) {
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

    if (x$fmodel[i] == 6) {
      if (x$fpm[i] == 1) {
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

    if (x$fmodel[i] == 7) {
      if (x$fpm[i] == 1) {
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
    if (x$fmodel[i] == 8) {
      if (x$fpm[i] == 1) {
        cat("Multiple extrema (IncMultExtreme) shaped", "\n")
        cat("Start increasing before omega_1 and alternate monotone shape at omega_j afterwards", "\n")
        cat("f(x) \t= int_0^x Z^2(s)h(s) ds - zeta + alpha_0 ", "\n")
        cat("zeta \t= min(0,min int_0^x Z^2(s)h(s)ds)", "\n")
      } else {
        cat("Multiple extrema (DecMultExtreme) shaped", "\n")
        cat("Start decreasing before omega_1 and alternate monotone shape at omega_j afterwards", "\n")
        cat("f(x) \t= -int_0^x Z^2(s)h(s) ds + zeta + alpha_0 ", "\n")
        cat("zeta \t= min(0,min int_0^x Z^2(s)h(s)ds)", "\n")
      }
      cat("h(x) \t= sum_{j=1}^{J} (-1)^(j-1) {1- exp[psi_j*(x-omega_j)]}/{exp[psi_j*(x-omega_j)]+1} ", "\n")
    }
    if (x$fmodel[i] > 2 & x$fmodel[i] < 7) {
      cat("Model for f includes linear term alpha", "\n")
    }
    cat("\n")
    cat("-----", "\n")
  }
  cat("Cosine basis ", "\n")
  cat("xrange   = xmax-xmin", "\n")
  cat("phi_0(x) = 1/sqrt(xmax-xmin) for xmin < x < xmax", "\n")
  cat("phi_k(x) = sqrt(2/(xmax-xmin))*cos(pi*k*(x-xmin)/(xmax-xmin))", "\n")
  cat("", "\n")
  if (x$model == "gbsar") {
    cat("Priors:", "\n")
    cat("beta    ~ N(b0,B0)", "\n")
    if (sum(x$fmodel > 2) > 0)
      cat("alpha   ~ N(m0,v0)I(delta*alpha>0)", "\n")
    if (sum(x$fmodel == 1) > 0)
      cat("theta_0 ~ N(0,v0)", "\n")
    if (sum(x$fmodel >= 2) > 0)
      cat("theta_0 ~ N(0,v0)I(theta_0 > 0)", "\n")
    if (sum(x$fmodel == 1) > 0)
      cat("theta_k ~ N(0,tau^2*exp(-gamma*k)", "\n")
    if (sum(x$fmodel >= 2) > 0)
      cat("theta_k ~ N(0,tau^2*exp(-gamma*k)", "\n")
  } else {
    cat("Scale invariant priors:", "\n")
    cat("beta|sigma \t~ N(b0,sigma2*B0)", "\n")
    if (sum(x$fmodel > 2) > 0)
      cat("alpha|sigma\t~ N(m0,sigma2*v0)I(delta*alpha>0)", "\n")
    if (sum(x$fmodel == 1) > 0)
      cat("theta_0|sigma \t~ N(0,sigma2*v0)", "\n")
    if (sum(x$fmodel >= 2) > 0)
      cat("theta_0|sigma\t~ N(0,sigma*v0)I(theta_0 > 0)", "\n")
    if (sum(x$fmodel == 1) > 0)
      cat("theta_k|sigma \t~ N(0,sigma2*tau^2*exp(-gamma*k)", "\n")
    if (sum(x$fmodel >= 2) > 0)
      cat("theta_k|sigma \t~ N(0,sigma*tau^2*exp(-gamma*k)", "\n")
  }
  cat("", "\n")
  cat("Smoothing parameters tau and gamma", "\n")
  cat("Choice of two priors for tau2", "\n")
  if (x$prior$iflagprior == 0)
    cat("T Smoother: tau2 ~ IG(r0/2,s0/2)", "\n")
  if (x$prior$iflagprior == 1)
    cat("Lasso Smoother: tau2 ~ Exp(u0)", "\n")
  cat("gamma ~ Exp(w0)", "\n")
  cat("\n")
  cat("Note: posterior of tau and gamma have banana contours", "\n")
  cat("zeta = ln(tau2) - kbar*gamma is less dependent with gamma or log(gamma)", "\n")
  cat("", "\n")
  if ((x$fmodel[i] == 5) || (x$fmodel[i] == 6) || (x$fmodel[i] == 7)) {
    cat("S models uses squish function (reparameterization of hyperbolic tangent)", "\n")
    cat("that depens on slope psi and location omega (inflection point for S).", "\n")
    cat("psi \t~ N(m0,v0)I(psi > 0)  Truncated Normal", "\n")
    cat("omega \t~ N(m0,v0)I(xmin < omega < xmax) and m0 = (xmin+xmax)/2", "\n")
    cat("\t", "\n")
  }
  if (x$fmodel[i] == 8) {
    cat("Multiple Extremes models uses squish functions (reparameterization of hyperbolic tangent)", "\n")
    cat("that depens on slope psi and location omega (inflection point for S).", "\n")
    cat("psi_j \t~ N(m0j,v0)I(psi_j > 0)  Truncated Normal", "\n")
    cat("omega_j \t~ N(m0j,v0)I(omega_{j-1} < omega_j < omega_{j+1}) and m0j = xmin + (xmax-xmin)*j/(J+2)", "\n")
    cat("\t", "\n")
  }
}
