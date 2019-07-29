\name{wage.union}
\alias{wage.union}
\docType{data}

\title{Wage-Union data}

\description{This data set contains 534 observations on 11 variables.}

\usage{data(wage.union)}

\format{
\describe{
\item{\code{education}}{number of years of education.}
\item{\code{south}}{indicator of living in southern region of U.S.A.}
\item{\code{sex}}{gender indicator: 0=male,1=female.}
\item{\code{experience}}{number of years of work experience.}
\item{\code{union}}{indicator of trade union membership: 0=non-member, 1=member.}
\item{\code{wage}}{wages in dollars per hour.}
\item{\code{age}}{age in years.}
\item{\code{race}}{1=black, 2=Hispanic, 3=white.}
\item{\code{occupation}}{1=management, 2=sales, 3=clerical, 4=service, 5=professional, 6=other.}
\item{\code{sector}}{0=other, 1=manufacturing, 2=construction.}
\item{\code{married}}{indicator of being married: 0=unmarried, 1=married.}
}
}

\references{

Berndt, E.R. (1991) \emph{The Practice of Econometrics}. New York: Addison-Wesley.

Ruppert, D., Wand, M.P. and Carroll, R.J. (2003) \emph{Semiparametric Regression}. Cambridge University Press.

}
\examples{
\dontrun{
	data(wage.union)
	pairs(wage.union)
}
}

\keyword{datasets}
