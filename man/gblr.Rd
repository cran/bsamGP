\name{gblr}
\alias{gblr}

\title{Generalized Bayesian Linear Models}
\description{
	This function fits a Bayesian generalized linear regression model.
}

\usage{

gblr(formula, data = NULL, family, link, mcmc = list(), prior = list(),
marginal.likelihood = TRUE, algorithm = c('AM', 'KS'), verbose = FALSE)
}

\arguments{
  \item{formula}{		an object of class \dQuote{\code{\link{formula}}} }

  \item{data}{	    an optional data frame.}

  \item{family}{		a description of the error distribution to be used in the model:
				The family contains bernoulli (\dQuote{bernoulli}),
				poisson (\dQuote{poisson}), negative-binomial (\dQuote{negative.binomial}),
				poisson-gamma mixture (\dQuote{poisson.gamma}). }

  \item{link}{			a description of the link function to be used in the model. }

  \item{mcmc}{			a list giving the MCMC parameters.
				The list includes the following integers (with default values in parentheses):
				\code{nblow (10000)} giving the number of MCMC in transition period,
				\code{nskip (10)} giving the thinning interval,
				\code{smcmc (1000)} giving the number of MCMC for analysis, and
				\code{ndisp (1000)} giving the number of saved draws to be displayed on screen
				(the function reports on the screen when every \code{ndisp} iterations have been carried out). }

  \item{prior}{			a list giving the prior information. The list includes the following parameters
				(default values specify the non-informative prior):
				\code{beta_m0} and \code{beta_v0} giving the hyperparameters of the
                     		multivariate normal distribution for parametric part including intercept,
                     		\code{kappa_m0} and \code{kappa_v0} giving the prior mean and variance of the
                     		gammal prior distribution for dispersion parameter (negative-binomial). }

  \item{marginal.likelihood}{	a logical variable indicating whether the log marginal likelihood is calculated.
				The methods of Gelfand and Dey (1994) is used. }

	\item{algorithm}{ a description of the algorithm to be used in the fitting of the logistic model:
	      The algorithm contains the Gibbs sampler based on the Kolmogorov-Smirnov distribution (\code{KS})
	      and an adaptive Metropolis  algorithm (\code{AM}).}

	 \item{verbose}{ a logical variable. If \code{TRUE}, the iteration number and the Metropolis acceptance rate are printed to the screen.}

}

\details{
  This generic function fits a Bayesian generalized linear regression models.

  Let \eqn{y_i} and \eqn{w_i} be the response and the vector of parametric predictors, respectively.
  The model is as follows.

  \deqn{y_i | \mu_i \sim F(\mu_i), }
  \deqn{g(\mu_i) = w_i^T\beta, ~ i=1,\ldots,n, }
  where \eqn{g(\cdot)} is a link function and \eqn{F(\cdot)} is a distribution of an exponential family.

  For unknown coefficients, the following prior is assumed for \eqn{\beta}:
  \deqn{\beta \sim N(m_{0,\beta}, V_{0,\beta})}

  The prior for the dispersion parameter of negative-binomial regression is
  \deqn{\kappa \sim Ga(r_0, s_0)}

}

\value{
  An object of class \code{blm} representing the generalized Bayesian linear model fit.
  Generic functions such as \code{print}, \code{fitted} and \code{plot} have methods to show the results of the fit.

  The MCMC samples of the parameters in the model are stored in the list \code{mcmc.draws},
  the posterior samples of the fitted values are stored in the list \code{fit.draws}, and
  the MCMC samples for the log marginal likelihood are saved in the list \code{loglik.draws}.
  The output list also includes the following objects:

  \item{post.est}{	posterior estimates for all parameters in the model. }

  \item{lmarg}{ 	log marginal likelihood using Gelfand-Dey method. }

  \item{family}{	the family object used. }

  \item{link}{		the link object used. }

  \item{methods}{	the method object used in the logit model. }

  \item{call}{		the matched call. }

  \item{mcmctime}{     running time of Markov chain from \code{system.time()}. }

}

\seealso{
\code{\link{blr}}, \code{\link{blq}}
}

\references{

Albert, J. H. and Chib, S. (1993) Bayesian Analysis of Binary and Polychotomous Response Data.
	\emph{Journal of the American Statistical Association}, \bold{88}, 669-679.

Holmes, C. C. and Held, L. (2006) Bayesian Auxiliary Variables Models for Binary and Multinomial Regression.
	\emph{Bayesian Analysis}, \bold{1}, 145-168.

Gelfand, A. E. and Dey, K. K. (1994) Bayesian Model Choice: Asymptotics and Exact Calculations.
	\emph{Journal of the Royal Statistical Society. Series B - Statistical Methodology}, \bold{56}, 501-514.

Roberts, G. O. and Rosenthal, J. S. (2009) Examples of Adaptive MCMC. \emph{Journal of Computational and Graphical Statistics}, \bold{18}, 349-367.

}

\examples{
############################
# Poisson Regression Model #
############################

# Simulate data
set.seed(1)

n <- 100
x <- runif(n)
y <- rpois(n, exp(0.5 + x*0.4))

# Fit the model with default priors and mcmc parameters
fout <- gblr(y ~ x, family = 'poisson', link = 'log')

# Summary
print(fout); summary(fout)

# Plot
plot(fout)

# fitted values
fitf <- fitted(fout)
}

\keyword{ generalized linear model }
\keyword{ parametric model }
