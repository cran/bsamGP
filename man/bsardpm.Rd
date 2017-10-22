\name{bsardpm}
\alias{bsardpm}

\title{Bayesian Shape-Restricted Spectral Analysis Regression with Dirichlet Process Mixture Errors}
\description{
	This function fits a Bayesian semiparametric regression model
	to estimate shape-restricted functions
	using a spectral analysis of Gaussian process priors.
	The model assumes that the errors follow a Dirichlet process mixture model.
}

\usage{

bsardpm(formula, xmin, xmax, nbasis, nint,
mcmc = list(), prior = list(), egrid, ngrid, location = TRUE,
shape = c('Free', 'Increasing', 'Decreasing', 'IncreasingConvex', 'DecreasingConcave',
'IncreasingConcave', 'DecreasingConvex', 'IncreasingS', 'DecreasingS',
'IncreasingRotatedS','DecreasingRotatedS','InvertedU','Ushape'))
}

\arguments{
  \item{formula}{		an object of class \dQuote{\code{\link{formula}}} }

  \item{xmin}{		a vector or scalar giving user-specific minimum values of x.
  				The default values are minimum values of x.}

  \item{xmax}{		a vector or scalar giving user-specific maximum values of x.
  				The default values are maximum values of x.}

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
				\code{theta0_m0} and \code{theta0_s0} giving the hyperparameters
				for prior distribution of the spectral coefficients
				(\code{theta0_m0} and \code{theta0_s0} are used when the functions have shape-restriction),
                     		\code{tau2_m0}, \code{tau2_s0} and \code{w0} giving the prior mean and standard deviation of
				smoothing prior (When iflagprior=1, tau2_m0 is only used as the hyperparameter),
				\code{beta_m0} and \code{beta_v0} giving the hyperparameters of the
                     		multivariate normal distribution for parametric part including intercept,
				\code{sigma2_m0} and \code{sigma2_v0} giving the prior mean and variance
				of the inverse gamma prior for the scale parameter of response,
                     		\code{alpha_m0} and \code{alpha_s0} giving the prior mean and standard deviation of the
                     		truncated normal prior distribution for the constant of integration,
				\code{iflagpsi} determining the prior of slope for logisitic function in \code{S} or \code{U} shaped
				(iflagpsi=1 (default), slope \eqn{\psi} is sampled and iflagpsi=0, \eqn{\psi} is fixed),
				\code{psifixed} giving initial value (iflagpsi=1) or fixed value (iflagpsi=0) of slope,
				\code{omega_m0} and \code{omega_s0} giving the prior mean and standard deviation of the
				truncated normal prior distribution for the inflection point of \code{S} or \code{U} shaped function.}

  \item{egrid}{			a vector giving grid points where the residual density estimate is evaluated.
  				The default range is from -10 to 10.}

  \item{ngrid}{			a vector giving number of grid points where the residual density estimate is evaluated.
  				The default value is 500.}

  \item{location}{		a logical value. If it is true, error density is modelled using location-scale mixture.}

  \item{shape}{			a vector giving types of shape restriction. }

}

\details{
  This generic function fits a Bayesian spectral analysis regression model
  for estimating shape-restricted functions using Gaussian process priors.
  For enforcing shape-restrictions, the model assumes that the derivatives of the functions are squares of Gaussian processes.
  The model also assumes that the errors follow a Dirichlet process mixture model.

  Let \eqn{y_i} and \eqn{w_i} be the response and the vector of parametric predictors, respectively.
  Further, let \eqn{x_{i,k}} be the covariate related to the response through an unknown shape-restricted function.
  The model for estimating shape-restricted functions is as follows.

  \deqn{y_i = w_i^T\beta + \sum_{k=1}^K f_k(x_{i,k}) + \epsilon_i, ~ i=1,\ldots,n, }
  where \eqn{f_k} is an unknown shape-restricted function of the scalar \eqn{x_{i,k} \in [0,1]} and
  the error terms \eqn{\{\epsilon_i\}} are a random sample from a Dirichlet process mixture model,

  1. scale mixture :
  \deqn{\epsilon_i \sim f(\epsilon) = \int N(\epsilon; 0,\sigma^2)dG(\sigma^2), }
  \deqn{G \sim DP(M,G0), ~~ G0 = Ga\left(\sigma^{-2}; \frac{r_{0,\sigma}}{2},\frac{s_{0,\sigma}}{2}\right).}

  2. location-scale mixture :
  \deqn{\epsilon_i \sim f(\epsilon) = \int N(\epsilon; \mu,\sigma^2)dG(\mu,\sigma^2), }
  \deqn{G \sim DP(M,G0), ~~ G0 = N\left(\mu;\mu_0,\kappa\sigma^2\right)Ga\left(\sigma^{-2}; \frac{r_{0,\sigma}}{2},\frac{s_{0,\sigma}}{2}\right).}

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

  For the spectral coefficients of functions without shape constraints, the scale-invariant prior is used
  (The intercept is included in \eqn{\beta}):
  \deqn{\theta_j | \tau, \gamma \sim N(0, \tau^2\exp[-j\gamma]), ~ j \ge 1}

  The priors for the spectral coefficients of shape restricted functions are:
  \deqn{\theta_0 \sim N(m_{\theta_0}, v^2_{\theta_0}), \quad
	\theta_j | \tau, \gamma \sim N(m_{\theta_j}, \tau^2\exp[-j\gamma]), ~ j \ge 1}

  To complete the model specification, the popular normal prior is assumed for \eqn{\beta}:
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

  \item{lpml}{		log pseudo marginal likelihood using Mukhopadhyay and Gelfand method. }

  \item{imodmet}{	the number of times to modify Metropolis. }

  \item{pmet}{		proportion of \eqn{\theta} accepted after burn-in. }

  \item{call}{		the matched call. }

  \item{mcmctime}{     running time of Markov chain from \code{system.time()}. }

}

\seealso{
\code{\link{bsar}}, \code{\link{bsaqdpm}}
}

\references{

Lenk, P. and Choi, T. (2017) Bayesian Analysis of Shape-Restricted Functions using
	Gaussian Process Priors. \emph{Statistica Sinica}, \bold{27}, 43-69.

MacEachern, S. N. and M{\"u}ller, P. (1998) Estimating mixture of Dirichlet process models.
	\emph{Journal of Computational and Graphical Statistics}, \bold{7}, 223-238.

Mukhopadhyay, S. and Gelfand, A. E. (1997) Dirichlet process mixed generalized linear models.
	\emph{Journal of the American Statistical Association}, \bold{92}, 633-639.

Neal, R. M. (2000) Markov chain sampling methods for Dirichlet process mixture models.
	\emph{Journal of Computational and Graphical Statistics}, \bold{9}, 249-265.



}

\examples{
\dontrun{
#####################
# Increasing-convex #
#####################

# Simulate data
set.seed(1)

n <- 200
x <- runif(n)
e <- c(rnorm(n/2, sd = 0.5), rnorm(n/2, sd = 3))
y <- exp(6*x - 3) + e

# Number of cosine basis functions
nbasis <- 50

# Fit the model with default priors and mcmc parameters
fout <- bsardpm(y ~ fs(x), nbasis = nbasis, shape = 'IncreasingConvex')

# Summary
print(fout); summary(fout)

# fitted values
fit <- fitted(fout)

# Plot
plot(fit, ask = TRUE)

}
}

\keyword{ Dirichlet processes }
\keyword{ Gaussian processes }
\keyword{ shape-restricted regression }
