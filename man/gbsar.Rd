\name{gbsar}
\alias{gbsar}

\title{Bayesian Shape-Restricted Spectral Analysis for Generalized Partial Linear Models}
\description{
	This function fits a Bayesian generalized partial linear regression model
	to estimate shape-restricted functions
	using a spectral analysis of Gaussian process priors.
}

\usage{

gbsar(formula, xmin, xmax, family, link, nbasis, nint, mcmc = list(), prior = list(),
shape = c('Free','Increasing','Decreasing','IncreasingConvex','DecreasingConcave',
          'IncreasingConcave','DecreasingConvex','IncreasingS','DecreasingS',
          'IncreasingRotatedS','DecreasingRotatedS','InvertedU','Ushape'),
marginal.likelihood = TRUE, algorithm = c('AM', 'KS'), verbose = FALSE)
}

\arguments{
  \item{formula}{		an object of class \dQuote{\code{\link{formula}}} }

  \item{xmin}{		a vector or scalar giving user-specific minimum values of x.
  				The default values are minimum values of x.}

  \item{xmax}{		a vector or scalar giving user-specific maximum values of x.
  				The default values are maximum values of x.}

  \item{family}{		a description of the error distribution to be used in the model:
				The family contains bernoulli (\dQuote{bernoulli}), poisson (\dQuote{poisson}),
				negative-binomial (\dQuote{negative.binomial}), poisson-gamma mixture (\dQuote{poisson.gamma}). }

  \item{link}{			a description of the link function to be used in the model. }

  \item{nbasis}{		number of cosine basis functions. }

  \item{nint}{			number of grid points where the unknown function is evaluated for plotting.
				The default is 200.}

  \item{mcmc}{			a list giving the MCMC parameters.
				The list includes the following integers (with default values in parentheses):
                     		\code{nblow0 (1000)} giving the number of initialization period for adaptive metropolis,
                     		\code{maxmodmet (5)} giving the maximum number of times to modify metropolis,
				\code{nblow (10000)} giving the number of MCMC in transition period,
				\code{nskip (10)} giving the thinning interval,
				\code{smcmc (1000)} giving the number of MCMC for analysis, and
                     		\code{ndisp (1000)} giving the number of saved draws to be displayed on screen
				(the function reports on the screen when every
				\code{ndisp} iterations have been carried out). }

  \item{prior}{			a list giving the prior information. The list includes the following parameters
				(default values specify the non-informative prior):
				\code{iflagprior} choosing a smoothing prior for spectral coefficients
				(iflagprior=0 assigns T-Smoother prior (default), iflagprior=1 chooses Lasso-Smoother prior),
				\code{theta_m0}, \code{theta0_m0} and \code{theta0_s0} giving the hyperparameters
				for prior distribution of the spectral coefficients
				(\code{theta0_m0} and \code{theta0_s0} are used when the functions have shape-restriction),
                     		\code{tau2_m0}, \code{tau2_s0} and \code{w0} giving the prior mean and standard deviation of
				smoothing prior (When iflagprior=1, tau2_m0 is only used as the hyperparameter),
				\code{beta_m0} and \code{beta_v0} giving the hyperparameters of the
                     		multivariate normal distribution for parametric part including intercept,
                     		\code{alpha_m0} and \code{alpha_s0} giving the prior mean and standard deviation of the
                     		truncated normal prior distribution for the constant of integration,
				\code{iflagpsi} determining the prior of slope for logisitic function in \code{S} or \code{U} shaped
				(iflagpsi=1 (default), slope \eqn{\psi} is sampled and iflagpsi=0, \eqn{\psi} is fixed),
				\code{psifixed} giving initial value (iflagpsi=1) or fixed value (iflagpsi=0) of slope,
				\code{omega_m0} and \code{omega_s0} giving the prior mean and standard deviation of the
				truncated normal prior distribution for the inflection point of \code{S} or \code{U} shaped function,
				\code{kappa_m0} and \code{kappa_v0} giving the prior mean and variance of the
                     		gammal prior distribution for dispersion parameter (negative-binomial). }

  \item{shape}{			a vector giving types of shape restriction. }

  \item{marginal.likelihood}{	a logical variable indicating whether the log marginal likelihood is calculated.
				The methods of Gelfand and Dey (1994) and Newton and Raftery (1994) are used.}

	\item{algorithm}{ a description of the algorithm to be used in the fitting of the logistic model:
	      The algorithm contains the Gibbs sampler based on the Kolmogorov-Smirnov distribution (\code{KS})
	      and an adaptive Metropolis algorithm (\code{AM}).}

	\item{verbose}{ a logical variable. If \code{TRUE}, the iteration number and the Metropolis acceptance rate are printed to the screen.}

}

\details{
  This generic function fits a Bayesian generalized partial linear regression models
  for estimating shape-restricted functions using Gaussian process priors.
  For enforcing shape-restrictions, they assumed that the derivatives of the functions are squares of Gaussian processes.

  Let \eqn{y_i} and \eqn{w_i} be the response and the vector of parametric predictors, respectively.
  Further, let \eqn{x_{i,k}} be the covariate related to the response through an unknown shape-restricted function.
  The model for estimating shape-restricted functions is as follows.

  \deqn{y_i | \mu_i \sim F(\mu_i), }
  \deqn{g(\mu_i) = w_i^T\beta + \sum_{k=1}^K f_k(x_{i,k}), ~ i=1,\ldots,n, }
  where \eqn{g(\cdot)} is a link function and \eqn{f_k} is an unknown nonlinear function of the scalar \eqn{x_{i,k} \in [0,1]}.

  The prior of function without shape restriction is:
  \deqn{f(x) = Z(x), }
  where \eqn{Z} is a second-order Gaussian process with mean function equal to zero and covariance function
  \eqn{\nu(s,t) = E[Z(s)Z(t)]} for \eqn{s, t \in [0, 1]}. The Gaussian process is expressed with
  the spectral representation based on cosine basis functions:
  \deqn{Z(x) = \sum_{j=0}^\infty \theta_j\varphi_j(x)}
  \deqn{\varphi_0(x) = 1 ~~ \code{and} ~~ \varphi_j(x) = \sqrt{2}\cos(\pi j x), ~ j \ge 1, ~ 0 \le x \le 1}

  The shape-restricted functions are modeled by assuming the \eqn{q}th derivatives of \eqn{f} are squares of Gaussian processes:
  \deqn{f^{(q)}(x) = \delta Z^2(x)h(x), ~~ \delta \in \{1, -1\}, ~~ q \in \{1, 2\},}
  where \eqn{h} is the squish function. For monotonic, monotonic convex, and concave functions, \eqn{h(x)=1}, while
  for \code{S} and \code{U} shaped functions, \eqn{h} is defined by
  \deqn{h(x) = \frac{1 - \exp[\psi(x - \omega)]}{1 + \exp[\psi(x - \omega)]}, ~~ \psi > 0, ~~ 0 < \omega < 1}

  For the spectral coefficients of functions without shape constraints, the following prior is used
  (The intercept is included in \eqn{\beta}):
  \deqn{\theta_j | \tau, \gamma \sim N(0, \tau^2\exp[-j\gamma]), ~ j \ge 1}

  The priors for the spectral coefficients of shape restricted functions are:
  \deqn{\theta_0 | \sim N(m_{\theta_0}, v^2_{\theta_0}), \quad
	\theta_j | \tau, \gamma \sim N(m_{\theta_j}, \tau^2\exp[-j\gamma]), ~ j \ge 1}

  To complete the model specification, the following prior is assumed for \eqn{\beta}:
  \deqn{\beta | \sim N(m_{0,\beta}, V_{0,\beta})}

}

\value{
  An object of class \code{bsam} representing the Bayesian spectral analysis model fit.
  Generic functions such as \code{print}, \code{fitted} and \code{plot} have methods to show the results of the fit.

  The MCMC samples of the parameters in the model are stored in the list \code{mcmc.draws},
  the posterior samples of the fitted values are stored in the list \code{fit.draws}, and
  the MCMC samples for the log marginal likelihood are saved in the list \code{loglik.draws}.
  The output list also includes the following objects:

  \item{post.est}{	posterior estimates for all parameters in the model. }

  \item{lmarg.gd}{	log marginal likelihood using Gelfand-Dey method. }

  \item{lmarg.nr}{	log marginal likelihood using Netwon-Raftery method, which is biased. }

  \item{family}{	the family object used. }

  \item{link}{		the link object used. }

  \item{call}{		the matched call. }

  \item{mcmctime}{      running time of Markov chain from \code{system.time()}. }

}

\seealso{
\code{\link{bsaq}}, \code{\link{bsar}}
}

\references{

Albert, J. H. and Chib, S. (1993) Bayesian Analysis of Binary and Polychotomous Response Data.
	\emph{Journal of the American Statistical Association}, \bold{88}, 669-679.

Holmes, C. C. and Held, L. (2006) Bayesian Auxiliary Variables Models for Binary and Multinomial Regression.
	\emph{Bayesian Analysis}, \bold{1}, 145-168.

Lenk, P. and Choi, T. (2017) Bayesian Analysis of Shape-Restricted Functions using
	Gaussian Process Priors. \emph{Statistica Sinica}, \bold{27}, 43-69.

Gelfand, A. E. and Dey, K. K. (1994) Bayesian model choice: asymptotics and exact calculations.
	\emph{Journal of the Royal Statistical Society. Series B - Statistical Methodology}, \bold{56}, 501-514.

Newton, M. A. and Raftery, A. E. (1994) Approximate Bayesian inference with the weighted likelihood bootstrap (with discussion).
	\emph{Journal of the Royal Statistical Society. Series B - Statistical Methodology}, \bold{56}, 3-48.

Roberts, G. O. and Rosenthal, J. S. (2009) Examples of Adaptive MCMC. \emph{Journal of Computational and Graphical Statistics}, \bold{18}, 349-367.

}

\examples{
\dontrun{
###########################
# Probit Regression Model #
###########################

# Simulate data
set.seed(1)

f <- function(x) 1.5 * sin(pi * x)

n <- 1000
b <- c(1,-1)
rho <- 0.7
u  <- runif(n, min = -1, max = 1)
x  <- runif(n, min = -1, max = 1)
w1 <- runif(n, min = -1, max = 1)
w2 <- round(f(rho * x + (1 - rho) * u))
w  <- cbind(w1, w2)

y  <- w \%*\% b + f(x) + rnorm(n)
y <- (y > 0)

# Number of cosine basis functions
nbasis <- 50

# Fit the model with default priors and mcmc parameters
fout <- gbsar(y ~ w1 + w2 + fs(x), family = "bernoulli", link = "probit",
              nbasis = nbasis, shape = 'Free')

# Summary
print(fout); summary(fout)

# fitted values
fit <- fitted(fout)

# Plot
plot(fit, ask = TRUE)

######################################
# Logistic Additive Regression Model #
######################################

# Wage-Union data
data(wage.union); attach(wage.union)

race[race==1 | race==2]=0
race[race==3]=1

y <- union
w <- cbind(race,sex,south)
x <- cbind(wage,education,age)

# mcmc parameters
mcmc <- list(nblow0 = 10000,
             nblow = 10000,
             nskip = 10,
             smcmc = 1000,
             ndisp = 1000,
             maxmodmet = 10)

foutGBSAR <- gbsar(y ~ race + sex + south + fs(wage) + fs(education) + fs(age),
                   family = 'bernoulli', link = 'logit', nbasis = 50, mcmc = mcmc,
                   shape = c('Free','Decreasing','Increasing'))

# fitted values
fitGBSAR <- fitted(foutGBSAR)

# Plot
plot(fitGBSAR, ask = TRUE)

}
\dontshow{
  # Toy Example :  Probit Regression Model #
# Simulate data
set.seed(1)

f <- function(x) 1.5 * sin(pi * x)

n <- 100
b <- c(1,-1)
rho <- 0.7
u  <- runif(n, min = -1, max = 1)
x  <- runif(n, min = -1, max = 1)
w1 <- runif(n, min = -1, max = 1)
w2 <- round(f(rho * x + (1 - rho) * u))
w  <- cbind(w1, w2)

y  <- w \%*\% b + f(x) + rnorm(n)
y <- (y > 0)

# Number of cosine basis functions
nbasis <- 10

# Fit the model with default priors and mcmc parameters
fout <- gbsar(y ~ w1 + w2 + fs(x), family = "bernoulli", link = "probit",
              nbasis = nbasis, shape = 'Free')

# Summary
print(fout); summary(fout)

# fitted values
fit <- fitted(fout)

# Plot
plot(fit, ask = FALSE)
}
}

\keyword{ Gaussian processes }
\keyword{ shape-restricted partial linear regression }
