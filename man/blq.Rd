\name{blq}
\alias{blq}

\title{Bayesian Quantile Regression}
\description{
	This function fits a Bayesian quantile regression model.
}

\usage{
blq(formula, data = NULL, p, mcmc = list(), prior = list(), marginal.likelihood = TRUE)
}

\arguments{
  \item{formula}{		an object of class \dQuote{\code{\link{formula}}} }

  \item{data}{	    an optional data frame.}

  \item{p}{			quantile of interest (default=0.5). }

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
				of the inverse gamma prior for the scale parameter of response. }

  \item{marginal.likelihood}{	a logical variable indicating whether the log marginal likelihood is calculated.
				The methods of Gelfand and Dey (1994) is used. }
}

\details{
  This generic function fits a Bayesian quantile regression model.

  Let \eqn{y_i} and \eqn{w_i} be the response and the vector of parametric predictors, respectively.
  Further, let \eqn{x_{i,k}} be the covariate related to the response, linearly.
  The model is as follows.

  \deqn{y_i = w_i^T\beta + \epsilon_i, ~ i=1,\ldots,n, }
  where the error terms \eqn{\{\epsilon_i\}} are a random sample from an asymmetric Laplace distribution, \eqn{ALD_p(0,\sigma^2)},
  which has the following probability density function:
  \deqn{ALD_p(\epsilon; \mu, \sigma^2) = \frac{p(1-p)}{\sigma^2}\exp\Big(-\frac{(x-\mu)[p - I(x \le \mu)]}{\sigma^2}\Big),}
  where \eqn{0 < p < 1} is the skew parameter, \eqn{\sigma^2 > 0} is the scale parameter, \eqn{-\infty < \mu < \infty} is
  the location parameter, and \eqn{I(\cdot)} is the indication function.

  The conjugate priors are assumed for \eqn{\beta} and \eqn{\sigma}:
  \deqn{\beta | \sigma \sim N(m_{0,\beta}, \sigma^2V_{0,\beta}), \quad \sigma^2 \sim IG\Big(\frac{r_{0,\sigma}}{2}, \frac{s_{0,\sigma}}{2}\Big)}

}

\value{
  An object of class \code{blm} representing the Bayesian parametric linear model fit.
  Generic functions such as \code{print} and \code{fitted} have methods to show the results of the fit.

  The MCMC samples of the parameters in the model are stored in the list \code{mcmc.draws},
  the posterior samples of the fitted values are stored in the list \code{fit.draws}, and
  the MCMC samples for the log marginal likelihood are saved in the list \code{loglik.draws}.
  The output list also includes the following objects:

  \item{post.est}{	posterior estimates for all parameters in the model. }

  \item{lmarg}{		log marginal likelihood using Gelfand-Dey method. }

  \item{rsquarey}{ 	correlation between \eqn{y} and \eqn{\hat{y}}. }

  \item{call}{		the matched call. }

  \item{mcmctime}{     running time of Markov chain from \code{system.time()}. }

}

\seealso{
\code{\link{blr}}, \code{\link{gblr}}
}

\references{


Gelfand, A. E. and Dey, K. K. (1994) Bayesian model choice: asymptotics and exact calculations.
	\emph{Journal of the Royal Statistical Society. Series B - Statistical Methodology}, \bold{56}, 501-514.

Kozumi, H. and Kobayashi, G. (2011) Gibbs sampling methods for Bayesian quantile regression.
	\emph{Journal of Statistical Computation and Simulation}, \bold{81}(11), 1565-1578.


}

\examples{
#####################
# Simulated example #
#####################

# Simulate data
set.seed(1)

n <- 100
w <- runif(n)
y <- 3 + 2*w + rald(n, scale = 0.8, p = 0.5)

# Fit median regression
fout <- blq(y ~ w, p = 0.5)

# Summary
print(fout); summary(fout)

# fitted values
fit <- fitted(fout)

# Plots
plot(fout)
}


\keyword{ quantile regression model }
\keyword{ parametric model }
