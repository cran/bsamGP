'print.bsad' <- function(x, ...)
{
    cat('\n')
	  cat('Call:\n')
    print(x$call)
    cat('\n')

	  kappaloop <- x$mcmc$kappaloop
	  nblow <- x$mcmc$nblow
	  nskip <- x$mcmc$nskip
	  smcmc <- x$mcmc$smcmc
	  nobs <- x$nobs
	  xmin <- x$xmin
	  xmax <- x$xmax
	  nint <- x$nint
	  MaxNCos <- x$MaxNCos

    cat( "Bayesian Spectral Analysis Density Estimation (BSAD)",'\n');
    cat( "- Density Estimation using Log-Gaussian Process -",'\n');
    cat('\n');
    if (x$parametric == 'normal'){
        cat( "Parametric model is normal",'\n');
    }else if (x$parametric == 'gamma'){
        cat( "Parametric model is gamma",'\n');
    }else if (x$parametric == 'laplace'){
    	cat( "Parametric model is laplace",'\n');
    }else{
    	cat( "Nonparametric model",'\n');
    }
    cat( "--------------------------------------------------------------------",'\n');
    cat('\n')
    cat( "# of cosines Kappa is random on 0, 1, ..., MaxNCos ",'\n');
    cat('\n')
    cat( "f(x) = exp[Z(x)] /int exp[Z(s)] ds",'\n');
    cat( "Discretize problem",'\n');
    cat( "f = exp(Y_i)/(sum exp(y_j) delta)",'\n');
    if(x$parametric == 'none'){
    	cat( "Use nonparametric regression model for Y",'\n');
    	cat( "Y = mu + Phi*Z*theta + epsilon with Z = 0 or 1",'\n');
    }else{
    	cat( "Use semiparametric regression model for Y",'\n');
    	cat( "Y = D*beta + Phi*Z*theta + epsilon with Z = 0 or 1",'\n');
    }
    cat( "[epsilon] = N(0,sigma^2I)",'\n');
    cat( "Smoothing prior on theta",'\n');
    cat( "Uses Slice Sampling to generate Y",'\n');
    cat('\n')
    cat( "--------------------------------------------------------------------",'\n');
    cat('\n')
    cat( "Number of observations   = ", x$nobs, sep='','\n');
    cat( "Estimation Interval: min = ", x$xmin, " and max = ", x$xmax, sep='','\n');
    cat( "Number of intervals      = ", x$nint, sep='','\n');
    cat( "Number of cosine terms   = ", x$MaxNCos, sep='','\n');
    cat('\n')
    cat( "--------------------------------------------------------------------",'\n');
    cat('\n')
    cat( "Number of transition iterations           = ", (1+kappaloop)*nblow, sep='','\n');
    cat( "Number of iterations for analysis         = ", smcmc, sep='','\n');
    cat( "Number of loops between generating kappa  = ", kappaloop, sep='','\n');
    cat( "Number of loops between saved interations = ", kappaloop*nskip, sep='','\n');
    cat( "Total number of MCMC iterations           = ", (1+kappaloop)*(nblow+smcmc), sep='','\n');
    cat('\n');
    if(x$parametric != 'none'){
    	cat( "********************************************************************",'\n');
    	cat('\n');
    	cat( "Posterior Probabilities of Parametric vs Semi-parametric with kappa",'\n')
    	bout = x$PostProbs
    	names(bout) <- c("Para",'SemiPara')
    	print(bout)
    }
    if(x$marginal.likelihood){
    	cat( "********************************************************************",'\n');
    	cat('\n');
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
    cat( "====================================================================",'\n');
    cat('\n');
    if(x$parametric!='none'){
    	cat( "Parametric Model",'\n');
    	cat('\n');
    	cat( "Regression Coefficients",'\n');
    	sout = c("PostMean","STD","2.5%","50%","97.5%" );
    	if (x$parametric=='normal'){
    		name <- c("X","X^2");
    	}else if (x$parametric=='gamma'){
        	name <- c("ln(X-Min)","X-Min");
    	}else if (x$parametric=='laplace'){
    		name <- c('|X|')
    	}
    	bout <- cbind(x$post.est$betaParm, x$post.est$betaPars)
    	bout <- cbind(bout, apply(x$mcmc.draws$betaPar, 2, function(x) quantile(x,probs = 0.025)),
    	              apply(x$mcmc.draws$betaPar, 2, function(x) quantile(x,probs = 0.5)),
    	              apply(x$mcmc.draws$betaPar, 2, function(x) quantile(x,probs = 0.975)))
    	colnames(bout) <- sout
    	rownames(bout) <- name
    	print(bout)
    	cat('\n');
    	cat( "====================================================================",'\n');
    	cat('\n');
    }
    cat( "Semiparametric Model: Estimation of Kappa",'\n');
    cat('\n');
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
    	bout=cbind(x$post.est$betam,x$post.est$betas)
    	bout=cbind(bout,apply(x$mcmc.draws$beta,2,function(x) quantile(x,probs=0.025)),
    					apply(x$mcmc.draws$beta,2,function(x) quantile(x,probs=0.5)),
    					apply(x$mcmc.draws$beta,2,function(x) quantile(x,probs=0.975)))
    	colnames(bout)=sout
    	rownames(bout)=name
    	print(bout)
    	cat('\n');
    	cat( "--------------------------------------------------------------------",'\n');
    	cat('\n');
    }
    cat( "Model Order kappa",'\n');
    sout = c("PostMean","STD","2.5%","50%","97.5%" );
    bout = c(mean(x$mcmc.draws$kappa),sd(x$mcmc.draws$kappa),
    		 quantile(x$mcmc.draws$kappa,probs=c(0.025,0.5,0.975)))
    names(bout)=sout
    print(bout)
    cat('\n');
    cat( "Smoothing Parameter tau",'\n');
    sout = c("PostMean","STD","2.5%","50%","97.5%" );
    bout = c(mean(x$mcmc.draws$tau),sd(x$mcmc.draws$tau),
    		 quantile(x$mcmc.draws$tau,probs=c(0.025,0.5,0.975)))
    names(bout)=sout
    print(bout)
    cat('\n');
    cat( "Smoothing parameter gamma",'\n');
    sout = c("PostMean","STD","2.5%","50%","97.5%" );
    bout = c(mean(x$mcmc.draws$gam),sd(x$mcmc.draws$gam),
    		 quantile(x$mcmc.draws$gam,probs=c(0.025,0.5,0.975)))
    names(bout)=sout
    print(bout)
    cat('\n');
    cat( "====================================================================",'\n');
    cat('\n');
    cat( "Semiparametric Model: Kappa = Max # of theta",'\n');
    cat('\n');
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
    	bout=cbind(x$post.est$betaMaxKappam,x$post.est$betaMaxKappas)
    	bout=cbind(bout,apply(x$mcmc.draws$betaMaxKappa,2,function(x) quantile(x,probs=0.025)),
    					apply(x$mcmc.draws$betaMaxKappa,2,function(x) quantile(x,probs=0.5)),
    					apply(x$mcmc.draws$betaMaxKappa,2,function(x) quantile(x,probs=0.975)))
    	colnames(bout)=sout
    	rownames(bout)=name
    	print(bout)
    	cat('\n');
    	cat( "--------------------------------------------------------------------",'\n');
    	cat('\n');
    }
    cat( "Smoothing Parameter tau using Max Kappa",'\n');
    sout = c("PostMean","STD","2.5%","50%","97.5%" );
    bout = c(mean(x$mcmc.draws$tauMaxKappa),sd(x$mcmc.draws$tauMaxKappa),
    		 quantile(x$mcmc.draws$tauMaxKappa,probs=c(0.025,0.5,0.975)))
    names(bout)=sout
    print(bout)
    cat('\n');
    cat( "Smoothing parameter gamma using Max Kappa",'\n');
    sout = c("PostMean","STD","2.5%","50%","97.5%" );
    bout = c(mean(x$mcmc.draws$gamMaxKappa),sd(x$mcmc.draws$gamMaxKappa),
    		 quantile(x$mcmc.draws$gamMaxKappa,probs=c(0.025,0.5,0.975)))
    names(bout)=sout
    print(bout)
    cat('\n');
    cat( "********************************************************************",'\n');
    cat('\n');
}
