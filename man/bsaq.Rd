\name{bsaq}
\alias{bsaq}

\title{Bayesian Shape-Restricted Spectral Analysis Quantile Regression}
\description{
	This function fits a Bayesian semiparametric quantile regression model
	to estimate shape-restricted functions
	using a spectral analysis of Gaussian process priors.
}

\usage{

bsaq(formula, xmin, xmax, p, nbasis, nint, mcmc = list(), prior = list(),
shape = c('Free', 'Increasing', 'Decreasing', 'IncreasingConvex', 'DecreasingConcave',
'IncreasingConcave', 'DecreasingConvex', 'IncreasingS', 'DecreasingS',
'IncreasingRotatedS', 'DecreasingRotatedS', 'InvertedU', 'Ushape'),
marginal.likelihood = TRUE, spm.adequacy = FALSE, verbose = FALSE)
}

\arguments{
  \item{formula}{		an object of class \dQuote{\code{\link{formula}}} }

  \item{xmin}{		a vector or scalar giving user-specific minimum values of x.
  				The default values are minimum values of x.}

  \item{xmax}{		a vector or scalar giving user-specific maximum values of x.
  				The default values are maximum values of x.}

  \item{p}{			quantile of interest (default=0.5). }

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

  \item{shape}{			a vector giving types of shape restriction. }

  \item{marginal.likelihood}{	a logical variable indicating whether the log marginal likelihood is calculated.
				The methods of Gelfand and Dey (1994) and Newton and Raftery (1994) are used. }

  \item{spm.adequacy}{		a logical variable indicating whether the log marginal likelihood of linear model is calculated.
				The marginal likelihood gives the values of the linear regression model excluding the nonlinear parts.}

	\item{verbose}{   a logical variable. If \code{TRUE}, the iteration number and the Metropolis acceptance rate are printed to the screen.}

}

\details{
  This generic function fits a Bayesian spectral analysis quantile regression model
  for estimating shape-restricted functions using Gaussian process priors.
  For enforcing shape-restrictions, the model assumed that the derivatives of the functions are squares of Gaussian processes.

  Let \eqn{y_i} and \eqn{w_i} be the response and the vector of parametric predictors, respectively.
  Further, let \eqn{x_{i,k}} be the covariate related to the response through an unknown shape-restricted function.
  The model for estimating shape-restricted functions is as follows.

  \deqn{y_i = w_i^T\beta + \sum_{k=1}^K f_k(x_{i,k}) + \epsilon_i, ~ i=1,\ldots,n, }
  where \eqn{f_k} is an unknown shape-restricted function of the scalar \eqn{x_{i,k} \in [0,1]} and
  the error terms \eqn{\{\epsilon_i\}} are a random sample from an asymmetric Laplace distribution, \eqn{ALD_p(0,\sigma^2)},
  which has the following probability density function:
  \deqn{ALD_p(\epsilon; \mu, \sigma^2) = \frac{p(1-p)}{\sigma^2}\exp\Big(-\frac{(x-\mu)[p - I(x \le \mu)]}{\sigma^2}\Big),}
  where \eqn{0 < p < 1} is the skew parameter, \eqn{\sigma^2 > 0} is the scale parameter, \eqn{-\infty < \mu < \infty} is
  the location parameter, and \eqn{I(\cdot)} is the indication function.

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
  \deqn{\theta_j | \sigma, \tau, \gamma \sim N(0, \sigma^2\tau^2\exp[-j\gamma]), ~ j \ge 1}

  The priors for the spectral coefficients of shape restricted functions are:
  \deqn{\theta_0 | \sigma \sim N(m_{\theta_0}, \sigma v^2_{\theta_0}), \quad
	\theta_j | \sigma, \tau, \gamma \sim N(m_{\theta_j}, \sigma\tau^2\exp[-j\gamma]), ~ j \ge 1}

  To complete the model specification, the conjugate priors are assumed for \eqn{\beta} and \eqn{\sigma}:
  \deqn{\beta | \sigma \sim N(m_{0,\beta}, \sigma^2V_{0,\beta}), \quad \sigma^2 \sim IG\left(\frac{r_{0,\sigma}}{2}, \frac{s_{0,\sigma}}{2}\right)}

}

\value{
  An object of class \code{bsam} representing the Bayesian spectral analysis model fit.
  Generic functions such as \code{print}, \code{fitted} and \code{plot} have methods to show the results of the fit.

  The MCMC samples of the parameters in the model are stored in the list \code{mcmc.draws},
  the posterior samples of the fitted values are stored in the list \code{fit.draws}, and
  the MCMC samples for the log marginal likelihood are saved in the list \code{loglik.draws}.
  The output list also includes the following objects:

  \item{post.est}{	posterior estimates for all parameters in the model. }

  \item{lmarg.lm}{	log marginal likelihood for linear quantile regression model. }

  \item{lmarg.gd}{	log marginal likelihood using Gelfand-Dey method. }

  \item{lmarg.nr}{	log marginal likelihood using Netwon-Raftery method, which is biased. }

  \item{rsquarey}{ 	correlation between \eqn{y} and \eqn{\hat{y}}. }

  \item{call}{		the matched call. }

  \item{mcmctime}{     running time of Markov chain from \code{system.time()}. }

}

\seealso{
\code{\link{bsar}}, \code{\link{gbsar}}
}

\references{

Lenk, P. and Choi, T. (2017) Bayesian Analysis of Shape-Restricted Functions using
	Gaussian Process Priors. \emph{Statistica Sinica}, 27: 43-69.

Gelfand, A. E. and Dey, K. K. (1994) Bayesian model choice: asymptotics and exact calculations.
	\emph{Journal of the Royal Statistical Society. Series B - Statistical Methodology}, \bold{56}, 501-514.

Kozumi, H. and Kobayashi, G. (2011) Gibbs sampling methods for Bayesian quantile regression.
	\emph{Journal of Statistical Computation and Simulation}, \bold{81}(11), 1565-1578.

Newton, M. A. and Raftery, A. E. (1994) Approximate Bayesian inference with the weighted likelihood bootstrap (with discussion).
	\emph{Journal of the Royal Statistical Society. Series B - Statistical Methodology}, \bold{56}, 3-48.


}

\examples{
\dontshow{
# Simulate data
set.seed(1)

n <- 100
x <- runif(n)
y <- log(1 + 10*x) + rald(n, scale = 0.5, p = 0.5)

# MCMC parameters
mcmc <- list(nblow = 1000, smcmc = 1000, nskip = 1)

# Fit the model with default priors and mcmc parameters
fout <- bsaq(y ~ fs(x), p = 0.5, nbasis = 30, mcmc = mcmc, shape = "Free")

# fitted values
fit <- fitted(fout)

# plots
plot(fit)
}

\dontrun{
######################
# Increasing-concave #
######################

# Simulate data
set.seed(1)

n <- 200
x <- runif(n)
y <- log(1 + 10*x) + rald(n, scale = 0.5, p = 0.5)

# Number of cosine basis functions
nbasis <- 50

# Fit the model with default priors and mcmc parameters
fout1 <- bsaq(y ~ fs(x), p = 0.25, nbasis = nbasis,
              shape = 'IncreasingConcave')
fout2 <- bsaq(y ~ fs(x), p = 0.5, nbasis = nbasis,
              shape = 'IncreasingConcave')
fout3 <- bsaq(y ~ fs(x), p = 0.75, nbasis = nbasis,
              shape = 'IncreasingConcave')

# fitted values
fit1 <- fitted(fout1)
fit2 <- fitted(fout2)
fit3 <- fitted(fout3)

# plots
plot(x, y, lwd = 2, xlab = 'x', ylab = 'y')
lines(fit1$xgrid, fit1$wbeta$mean[1] + fit1$fxgrid$mean, lwd=2, col=2)
lines(fit2$xgrid, fit2$wbeta$mean[1] + fit2$fxgrid$mean, lwd=2, col=3)
lines(fit3$xgrid, fit3$wbeta$mean[1] + fit3$fxgrid$mean, lwd=2, col=4)
legend('topleft', legend = c('1st Quartile', '2nd Quartile', '3rd Quartile'),
       lwd = 2, col = 2:4, lty = 1)

}
}

\keyword{ Gaussian processes }
\keyword{ shape-restricted quantile regression }
