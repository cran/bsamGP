\name{fitted.bsad}
\alias{fitted.bsad}
\title{ Compute fitted values for a bsad object}
\description{
  Computes pointwise posterior means and \eqn{100(1-\alpha)}\% credible intervals of the fitted Bayesian spectral analysis density estimation model.
}
\usage{
%fitted.bsad(object, alpha = 0.05, HPD = TRUE, \dots)
\method{fitted}{bsad}(object, alpha = 0.05, HPD = TRUE, \dots)
}
\arguments{
  \item{object}{ a \code{bsad} object}
  
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
A list object of class \code{fitted.bsad} containing posterior means and \eqn{100(1-\alpha)}\% credible intervals.
Generic function \code{plot} displays the results of the fit.

The output list includes the following objects:

\item{fpar}{	posterior estimates for parametric model.}

\item{fsemi}{	posterior estimates for semiparametric model.}

\item{fsemiMaxKappa}{	posterior estimates for semiparametric model with maximum number of basis.}

}

\seealso{
\code{\link{bsad}}
}

\examples{
## See examples for bsad
}

\keyword{ logistic Gaussian process}
\keyword{ density estimation }
