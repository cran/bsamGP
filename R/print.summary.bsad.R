'print.summary.bsad' <- function(x, ...)
{
  cat('\n')
  cat( "Number of observations   = ", x$nobs, sep='','\n');
  cat( "Estimation Interval: min = ", x$xmin, " and max = ", x$xmax, sep='','\n');
  cat( "Number of intervals      = ", x$nint, sep='','\n');
  cat( "Number of cosine terms   = ", x$MaxNCos, sep='','\n');
  cat('\n')
  cat( "Number of transition iterations           = ", (1+x$kappaloop)*x$nblow, sep='','\n');
  cat( "Number of iterations for analysis         = ", x$smcmc, sep='','\n');
  cat( "Number of loops between generating kappa  = ", x$kappaloop, sep='','\n');
  cat( "Number of loops between saved interations = ", x$kappaloop*x$nskip, sep='','\n');
  cat( "Total number of MCMC iterations           = ", (1+x$kappaloop)*(x$nblow+x$smcmc), sep='','\n');
  cat('\n');
  if(x$parametric != 'none'){
    cat( "Posterior Probabilities of Parametric vs Semi-parametric with kappa",'\n')
    bout = x$PostProbs
    names(bout) <- c("Para",'SemiPara')
    print(bout)
    cat('\n')
  }
  if(x$marginal.likelihood){
    cat( "Ln Marginal Distribution of the Data",'\n');
    if(x$parametric != 'none'){
      name <- c("Para","Semi","SemiMK");
      bout <- x$lmarg
      names(bout)=name
    }else{
      name <- c("Semi","SemiMK");
      bout <- x$lmarg
      names(bout)=name
    }
    print(bout)
    cat('\n');
  }
  if(x$parametric!='none'){
    cat("-----",'\n')
    cat( paste("Parametric Model : ", "\'", x$parametric, "\'", sep = ''),'\n');
    cat( "Regression Coefficients",'\n');
    sout = c("PostMean","STD","2.5%","50%","97.5%" );
    if (x$parametric=='normal'){
      name <- c("X","X^2");
    }else if (x$parametric=='gamma'){
      name <- c("ln(X-Min)","X-Min");
    }else if (x$parametric=='laplace'){
      name <- c('|X|')
    }
    bout <- cbind(x$betaParm, x$betaPars)
    bout <- cbind(bout, x$betaParql, x$betaParq2, x$betaParqu)
    colnames(bout) <- sout
    rownames(bout) <- name
    print(bout)
    cat('\n');
  }
  cat("-----",'\n')
  cat( "Semiparametric Model: Estimation of Kappa",'\n');
  if(x$parametric != 'none'){
    cat( "Regression Coefficients",'\n');
    sout <- c("PostMean","STD","2.5%","50%","97.5%" );
    if (x$parametric == 'normal'){
      name <- c("X","X^2");
    }else if (x$parametric == 'gamma'){
      name <- c("ln(X-Min)","X-Min");
    }else if (x$parametric == 'laplace'){
      name <- c('|X|')
    }
    bout=cbind(x$betam,x$betas)
    bout=cbind(bout,x$betaql, x$betaq2, x$betaqu)
    colnames(bout)=sout
    rownames(bout)=name
    print(bout)
    cat('\n');
  }
  cat( "Model Order kappa",'\n');
  sout = c("PostMean","STD","2.5%","50%","97.5%" );
  bout = c(x$kappam, x$kappas, x$kappaq)
  names(bout)=sout
  print(bout)
  cat('\n');
  cat( "Smoothing Parameter tau",'\n');
  sout = c("PostMean","STD","2.5%","50%","97.5%" );
  bout = c(x$taum, x$taus, x$tauq)
  names(bout)=sout
  print(bout)
  cat('\n');
  cat( "Smoothing parameter gamma",'\n');
  sout = c("PostMean","STD","2.5%","50%","97.5%" );
  bout = c(x$gammam, x$gammas, x$gammaq)
  names(bout)=sout
  print(bout)
  cat('\n');
  cat("-----",'\n')
  cat( "Semiparametric Model: Kappa = Max # of theta",'\n');
  if(x$parametric!='none'){
    cat( "Regression Coefficients",'\n');
    sout = c("PostMean","STD","2.5%","50%","97.5%" );
    if (x$parametric=='normal'){
      name = c("X","X^2");
    }else if (x$parametric=='gamma'){
      name = c("ln(X-Min)","X-Min");
    }else if (x$parametric=='laplace'){
      name = c('|X|')
    }
    bout=cbind(x$betaMaxKappam,x$betaMaxKappas)
    bout=cbind(bout,x$betaMaxKappaql, x$betaMaxKappaq2, x$betaMaxKappaqu)
    colnames(bout)=sout
    rownames(bout)=name
    print(bout)
    cat('\n');
  }
  cat( "Smoothing Parameter tau using Max Kappa",'\n');
  sout = c("PostMean","STD","2.5%","50%","97.5%" );
  bout = c(x$tauMaxKappam, x$tauMaxKappas, x$tauMaxKappaq)
  names(bout)=sout
  print(bout)
  cat('\n');
  cat( "Smoothing parameter gamma using Max Kappa",'\n');
  sout = c("PostMean","STD","2.5%","50%","97.5%" );
  bout = c(x$gamMaxKappam, x$gamMaxKappas, x$gamMaxKappaq)
  names(bout)=sout
  print(bout)
  cat('\n');
}
