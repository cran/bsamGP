\name{predict.bsam}
\alias{predict.bsam}
\title{Predict method for a bsam object}
\description{
  Computes the predicted values of Bayesian spectral analysis models.
}
\usage{
\method{predict}{bsam}(object, newp, newnp, alpha = 0.05, HPD = TRUE, type = "response", \dots)
}
\arguments{
  \item{object}{ a \code{bsam} object}

  \item{newp}{   an optional data of parametric components with which to predict. If omitted, the fitted values are returned.}

  \item{newnp}{  an optional data of nonparametric components with which to predict. If omitted, the fitted values are returned.}

  \item{alpha}{	 a numeric scalar in the interval (0,1) giving the \eqn{100(1-\alpha)}\% credible intervals.}

  \item{HPD}{    a logical variable indicating whether the \eqn{100(1-\alpha)}\% Highest Posterior Density (HPD) intervals are calculated.
		  If \code{HPD}=FALSE, the \eqn{100(1-\alpha)}\% equal-tail credible intervals are calculated.
		  The default is TRUE.}

  \item{type}{  the type of prediction required. \code{type = "response"} gives the posterior predictive samples as default. The \code{"mean"} option returns expectation of the posterior estimates. }

  \item{\dots}{  not used }
}
\details{
None.
}
\value{
A list object of class \code{predict.bsam} containing posterior means and \eqn{100(1-\alpha)}\% credible intervals.

The output list includes the following objects:

\item{fxobs}{	posterior estimates for unknown functions over observation.}

\item{wbeta}{	posterior estimates for parametric part.}

\item{yhat}{	posterior estimates for fitted values of either response or expectation of response.
		For \code{\link{gbsar}}, it gives posterior estimates for expectation of response.}

\item{fxResid}{ posterior estimates for fitted parametric residuals. Not applicable for \code{\link{gbsar}}.}

}

\seealso{
\code{\link{bsaq}}, \code{\link{bsar}}, \code{\link{gbsar}}
}

\examples{
\dontrun{

##########################################
# Increasing Convex to Concave (S-shape) #
##########################################

# simulate data
f <- function(x) 5*exp(-10*(x - 1)^4) + 5*x^2

set.seed(1)

n <- 100
x <- runif(n)
y <- f(x) + rnorm(n, sd = 1)

# Number of cosine basis functions
nbasis <- 50

# Fit the model with default priors and mcmc parameters
fout <- bsar(y ~ fs(x), nbasis = nbasis, shape = 'IncreasingConvex',
             spm.adequacy = TRUE)

# Prediction
xnew <- runif(n)
predict(fout, newnp = xnew)
}
}

