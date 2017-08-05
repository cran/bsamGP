\name{fitted.blm}
\alias{fitted.blm}
\title{ Compute fitted values for a blm object}
\description{
  Computes pointwise posterior means and 95\% credible intervals of the fitted Bayesian linear models.
}
\usage{
%fitted.blm(object, alpha = 0.05, HPD = TRUE, \dots)
\method{fitted}{blm}(object, alpha = 0.05, HPD = TRUE, \dots)
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
A list containing posterior means and 95\% credible intervals.

The output list includes the following object:

\item{wbeta}{	posterior estimates for regression function.}

}

\seealso{
\code{\link{blq}}, \code{\link{blr}}
}

\references{

Chen, M., Shao, Q. and Ibrahim, J. (2000) \emph{Monte Carlo Methods in Bayesian computation}. 
	Springer-Verlag New York, Inc.

}

\examples{
## See examples for blq and blr
}

