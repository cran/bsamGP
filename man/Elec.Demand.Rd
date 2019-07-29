\name{Elec.demand}
\alias{Elec.demand}
\docType{data}
\title{Electricity demand data}
\description{
  The Elec.demand data consists of 288 quarterly observations in Ontario from 1971 to 1994.
}
\usage{data(Elec.demand)}
\format{
  A data frame with 288 observations on the following 7 variables.
  \describe{
    \item{quarter}{date (yyyy-mm) from 1971 to 1994}
    \item{enerm}{electricity demand.}
    \item{gdp}{gross domestic product.}
    \item{pelec}{price of electricity.}
    \item{pgas}{price of natural gas.}
    \item{hddqm}{the number of heating degree days relative to a reference temperature.}
    \item{cddqm}{the number of cooling degree days relative to a reference temperature.}
  }
}
\source{
  Yatchew, A. (2003). \emph{Semiparametric Regression for the Applied Econometrician}.
	Cambridge University Press.
}
\references{
  Engle, R. F., Granger, C. W. J., Rice, J. and Weiss, A. (1986). Semiparametric estimates of
	the relation between weather and electricity sales.
	\emph{Journal of the American Statistical Association}, \bold{81}, 310-320.

  Lenk, P. and Choi, T. (2017). Bayesian analysis of shape-restricted functions using
	Gaussian process priors. \emph{Statistica Sinica}, \bold{27}, 43-69.

}
\examples{
\dontrun{
	data(Elec.demand)
	plot(Elec.demand)
}
}
\keyword{datasets}
