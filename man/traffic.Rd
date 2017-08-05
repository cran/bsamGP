\name{traffic}
\alias{traffic}
\docType{data}

\title{Monthly traffic accidents data}

\description{This data set contains 108 observations on 6 variables.}

\usage{data(traffic)}

\format{
\describe{
\item{\code{ln_number}}{logarithm of the number of monthly automobile accidents in the state of Michigan.}
\item{\code{month}}{months from January 1st, 1979 to Decembe 31st, 1987.}
\item{\code{ln_unemp}}{logarithm of unemployment rate}
\item{\code{spring}}{indicator for spring season.}
\item{\code{summer}}{indicator for summer season.}
\item{\code{fall}}{indicator for fall season.}
}
}

\references{

	Lenk (1999) Bayesian inference for semiparametric regression using a Fourier representation. 
	\emph{Journal of the Royal Statistical Society: Series B}, \bold{61}(4), 863-879. 

}
\examples{
\dontrun{
	data(traffic)
	pairs(traffic)
}
}

\keyword{datasets}
