\name{predict.blm}
\alias{predict.blm}
\title{ Predict method for a blm object}
\description{
  Computes predicted values of Bayesian linear models.
}
\usage{
	\method{predict}{blm}(object, newdata, alpha = 0.05, HPD = TRUE, \dots)
}
\arguments{
  \item{object}{  a \code{bsam} object}

  \item{newdata}{ an optional data matrix or vector with which to predict. If omitted, the fitted values are returned.}

  \item{alpha}{	  a numeric scalar in the interval (0,1) giving the \eqn{100(1-\alpha)}\% credible intervals.}

  \item{HPD}{     a logical variable indicating whether the \eqn{100(1-\alpha)}\% Highest Posterior Density (HPD) intervals are calculated.
		 If \code{HPD}=FALSE, the \eqn{100(1-\alpha)}\% equal-tail credible intervals are calculated.
		 The default is TRUE.}

  \item{\dots}{  not used }
}
\details{
None.
}
\value{
A list containing posterior means and 95\% credible intervals.

The output list includes the following objects:

\item{wbeta}{	posterior estimates for regression function.}

\item{yhat}{	posterior estimates for generalised regression function.}

}

\seealso{
\code{\link{blq}}, \code{\link{blr}}, \code{\link{gblr}}
}

\references{

Chen, M., Shao, Q. and Ibrahim, J. (2000) \emph{Monte Carlo Methods in Bayesian computation}.
	Springer-Verlag New York, Inc.

}

\examples{
\dontrun{
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

	  # Predict
	  new <- rnorm(n)
	  predict(fout, newdata = new)
}
}

