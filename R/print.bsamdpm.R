"print.bsamdpm" <- function(x, ...) {
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")

  if (x$model == "bsaqdpm") {
    cat("Bayesian Spectral Analysis Quantile Regression with Dirichlet process mixture errors \n")
    cat("\n")
    cat("Model:", "\n")
    cat("Y = w'beta + sum_k f_k(x_k) + int ALD(p;0,sigma^2)dG(sigma^2)", "\n")
    cat("G ~ DP(M,G0), G0=Ga(1/sigma^2; r0/2,s0/2)", "\n")
  } else {
    cat("Bayesian Spectral Analysis Regression with Dirichlet process mixture errors \n")
    cat("\n")
    cat("Model:", "\n")
    if (x$location) {
      cat("Y = w'beta + sum_k f_k(x_k) + int N(mu,sigma^2)dG(mu,sigma^2)", "\n")
      cat("G ~ DP(M,G0),", "\n")
      cat("G0 = N(mu;0,kappa/sigma^2)Ga(1/sigma^2; r0/2,s0/2)", "\n")
    } else {
      cat("Y = w'beta + sum_k f_k(x_k) + int N(0,sigma^2)dG(sigma^2)", "\n")
      cat("G ~ DP(M,G0), G0 = Ga(1/sigma^2; r0/2,s0/2)", "\n")
    }
  }
  cat('\n')
  cat("xmin < x < xmax", "\n")
  cat("Normalize f: int_a^b f(x)dx = 0", "\n")
  cat("beta has the Intercept", "\n")
  cat("Define", "\n")
  cat("Z(x) = sum_j=0^J theta_j*phi_j(x)", "\n")
  cat("g^a(x) = int_0^x |Z(s)|^2 dx", "\n")
  cat("       = theta' Phi^a(x)*theta", "\n")
  cat("where Phi^a(x) is matrix with ", "\n")
  cat("phi^a_{j,k}(x) = int_0^x phi_j(s) phi_k(s) ds in (j,k) ", "\n")
  cat("g^b(x) = int_0^x g^a(s) ds", "\n")
  cat("       = theta' Phi^b(x)*theta", "\n")
  cat("where Phi^b(x) is a matrix with", "\n")
  cat("phi^b_{j,k}(x) = int_0^x phi^a_{j,k}(s) ds", "\n")
  cat("\n")
}
