\name{blr}
\alias{blr}

\title{Bayesian Linear Regression}
\description{
	This function fits a Bayesian linear regression model using scale invariant prior.
}

\usage{
blr(formula, data = NULL, mcmc = list(), prior = list(), marginal.likelihood = TRUE)
}

\arguments{
  \item{formula}{		an object of class \dQuote{\code{\link{formula}}} }

  \item{data}{	    an optional data frame.}

  \item{mcmc}{			a list giving the MCMC parameters.
				The list includes the following integers (with default values in parentheses):
				\code{nblow (1000)} giving the number of MCMC in transition period,
				\code{nskip (1)} giving the thinning interval,
				\code{smcmc (1000)} giving the number of MCMC for analysis. }

  \item{prior}{			a list giving the prior information. The list includes the following parameters
				(default values specify the non-informative prior):
				\code{beta_m0} and \code{beta_v0} giving the hyperparameters of the
                     		multivariate normal distribution for parametric part including intercept,
				\code{sigma2_m0} and \code{sigma2_v0} giving the prior mean and variance
				of the inverse gamma prior for the scale parameter of response.}

  \item{marginal.likelihood}{	a logical variable indicating whether the log marginal likelihood is calculated. }

}

\details{
  This generic function fits a Bayesian linear regression model using scale invariant prior.

  Let \eqn{y_i} and \eqn{w_i} be the response and the vector of parametric predictors, respectively.
  The model for regression function is as follows.

  \deqn{y_i = w_i^T\beta + \epsilon_i, ~ i=1,\ldots,n, }
  where the error terms \eqn{\{\epsilon_i\}} are a random sample from a normal distribution, \eqn{N(0,\sigma^2)}.

  The conjugate priors are assumed for \eqn{\beta} and \eqn{\sigma}:
  \deqn{\beta | \sigma \sim N(m_{0,\beta}, \sigma^2V_{0,\beta}), \quad \sigma^2 \sim IG\Big(\frac{r_{0,\sigma}}{2}, \frac{s_{0,\sigma}}{2}\Big)}

}

\value{
  An object of class \code{blm} representing the Bayesian spectral analysis model fit.
  Generic functions such as \code{print} and \code{fitted} have methods to show the results of the fit.

  The MCMC samples of the parameters in the model are stored in the list \code{mcmc.draws} and
  the posterior samples of the fitted values are stored in the list \code{fit.draws}.
  The output list also includes the following objects:

  \item{post.est}{	posterior estimates for all parameters in the model. }

  \item{lmarg}{		log marginal likelihood. }

  \item{rsquarey}{ 	correlation between \eqn{y} and \eqn{\hat{y}}. }

  \item{call}{		the matched call. }

  \item{mcmctime}{      running time of Markov chain from \code{system.time()}. }

}

\seealso{
\code{\link{blq}}, \code{\link{gblr}}
}


\examples{
#####################
# Simulated example #
#####################

# Simulate data
set.seed(1)

n <- 100
w <- runif(n)
y <- 3 + 2*w + rnorm(n, sd = 0.8)

# Fit the model with default priors and mcmc parameters
fout <- blr(y ~ w)

# Summary
print(fout); summary(fout)

# Fitted values
fit <- fitted(fout)

# Plots
plot(fout)
}

\keyword{ linear regression model }
\keyword{ parametric model }
