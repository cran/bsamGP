\name{predict.bsamdpm}
\alias{predict.bsamdpm}
\title{Predict method for a bsamdpm object}
\description{
  Computes the predicted values of Bayesian spectral analysis models with Dirichlet process mixture errors.
}
\usage{
\method{predict}{bsamdpm}(object, newp, newnp, alpha = 0.05, HPD = TRUE, \dots)
}
\arguments{
  \item{object}{ a \code{bsamdpm} object}

  \item{newp}{   an optional data of parametric components with which to predict. If omitted, the fitted values are returned.}

  \item{newnp}{  an optional data of nonparametric components with which to predict. If omitted, the fitted values are returned.}

  \item{alpha}{	 a numeric scalar in the interval (0,1) giving the \eqn{100(1-\alpha)}\% credible intervals.}

  \item{HPD}{    a logical variable indicating whether the \eqn{100(1-\alpha)}\% Highest Posterior Density (HPD) intervals are calculated.
		  If \code{HPD}=FALSE, the \eqn{100(1-\alpha)}\% equal-tail credible intervals are calculated.
		  The default is TRUE.}

  \item{\dots}{  not used }
}
\details{
None.
}
\value{
A list object of class \code{predict.bsamdpm} containing posterior means and \eqn{100(1-\alpha)}\% credible intervals.

The output list includes the following objects:

\item{fxobs}{	posterior estimates for unknown functions over observation.}

\item{wbeta}{	posterior estimates for parametric part.}

\item{yhat}{	posterior estimates for fitted values of response.}

}

\seealso{
\code{\link{bsaqdpm}}, \code{\link{bsardpm}}
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

# Prediction
xnew <- runif(n)
predict(fout, newnp = xnew)

}
}

