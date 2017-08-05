\name{fitted.bsam}
\alias{fitted.bsam}
\title{ Compute fitted values for a bsam object}
\description{
  Computes pointwise posterior means and \eqn{100(1-\alpha)}\% credible intervals of the fitted Bayesian spectral analysis models.
}
\usage{
%fitted.bsam(object, alpha = 0.05, HPD = TRUE, \dots)
\method{fitted}{bsam}(object, alpha = 0.05, HPD = TRUE, \dots)
}
\arguments{
  \item{object}{ a \code{bsam} object}

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
A list object of class \code{fitted.bsam} containing posterior means and \eqn{100(1-\alpha)}\% credible intervals.
Generic function \code{plot} displays the results of the fit.

The output list includes the following objects:

\item{fxobs}{	posterior estimates for unknown functions over observation.}

\item{fxgrid}{	posterior estimates for unknown functions over grid points.}

\item{wbeta}{	posterior estimates for parametric part.}

\item{yhat}{	posterior estimates for fitted values of response.
		For \code{\link{gbsar}}, it gives posterior estimates for expectation of response.}

}

\seealso{
\code{\link{bsaq}}, \code{\link{bsaqdpm}}, \code{\link{bsar}}, \code{\link{bsardpm}}
}

\examples{
## See examples for bsaq, bsaqdpm, bsar, and bsardpm
}

