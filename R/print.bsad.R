'print.bsad' <- function(x, ...)
{
    cat('\n')
	  cat('Call:\n')
    print(x$call)
    cat('\n')

    cat("Bayesian Spectral Analysis Density Estimation (BSAD)",'\n');
    cat("- Density Estimation using Log-Gaussian Process -",'\n');
    cat('\n');
    cat('Model:')
    cat('\n')
    cat("f(x) = exp[Z(x)] /int exp[Z(s)] ds",'\n');
    cat("\n")
    cat("1. Discretize problem",'\n');
    cat("f = exp(Y_i)/(sum exp(Y_j) delta)",'\n');
    cat("2. Uses Slice Sampling to generate Y",'\n');
    if (x$parametric == 'none'){
      cat("3. Use nonparametric regression model for Y",'\n');
      cat("Y = mu + Phi*Z*theta + epsilon with Z = 0 or 1",'\n');
    }else{
      cat("3. Use semiparametric regression model for Y",'\n');
      cat("Y = D*beta + Phi*Z*theta + epsilon with Z = 0 or 1",'\n');
    }
    cat("4. [epsilon] = N(0,sigma^2I)",'\n');
    cat('\n')
    cat("Priors:",'\n')
    cat("theta_k|sigma 	~ N(0,tau^2*exp(-gamma*c_k)",'\n')
    cat('\n')
    cat('Smoothing parameters tau and gamma','\n')
    cat('tau2 ~ IG(r0/2, s0/2)','\n')
    cat('gamma ~ IG(u0/2, v0/2)','\n')
    cat('\n')
    cat('Choice of two possible parametrizations of c_k','\n')
    if (x$smooth == "geometric") {
      cat('Geometric: c_k = k', '\n')
    } else {
      cat('Algebraic: c_k = log(k + 1)','\n')
    }
    cat('\n')
}
