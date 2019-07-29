\name{fitted.bsamdpm}
\alias{fitted.bsamdpm}
\title{ Compute fitted values for a bsamdpm object}
\description{
  Computes pointwise posterior means and \eqn{100(1-\alpha)}\% credible intervals of the fitted Bayesian spectral analysis models with Dirichlet process mixture error.
}
\usage{
%fitted.bsamdpm(object, alpha = 0.05, HPD = TRUE, \dots)
\method{fitted}{bsamdpm}(object, alpha = 0.05, HPD = TRUE, \dots)
}
\arguments{
  \item{object}{ a \code{bsamdpm} object}

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
A list object of class \code{fitted.bsamdpm} containing posterior means and 95\% credible intervals.
Generic function \code{plot} displays the results of the fit.

The output list includes the following objects:

\item{edens}{	posterior estimate for unknown error distribution over grid points.}

\item{fxobs}{	posterior estimates for unknown functions over observation.}

\item{fxgrid}{	posterior estimates for unknown functions over grid points.}

\item{wbeta}{	posterior estimates for parametric part.}

\item{yhat}{	posterior estimates for fitted values of response.}

}

\seealso{
\code{\link{bsaqdpm}}, \code{\link{bsardpm}}
}

\examples{
## See examples for bsaqdpm and bsardpm
}

