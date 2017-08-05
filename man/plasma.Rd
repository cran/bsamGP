\name{plasma}
\alias{plasma}
\docType{data}

\title{A Data Set for Plasma Levels of Retinol and Beta-Carotene}

\description{This data set contains 314 observations on 14 variables.}

\usage{data(plasma)}

\format{
\describe{
\item{\code{age}}{Age (years).}
\item{\code{sex}}{Sex (1=Male, 2=Female).}
\item{\code{smoke}}{Smoking status (1=Never, 2=Former, 3=Current Smoker).}
\item{\code{vmi}}{BMI values (weight/(height^2)).}
\item{\code{vitas}}{Vitamin use (1=Yes,fairly often, 2=Yes, not often, 3=No).}
\item{\code{calories}}{Number of calories consumed per day.}
\item{\code{fat}}{Grams of fat consumed per day.}
\item{\code{fiber}}{Grams of fiber consumed per day.}
\item{\code{alcohol}}{Number of alcoholic drinks consumed per week.}
\item{\code{cholesterol}}{Cholesterol consumed (mg per day).}
\item{\code{beta diet}}{Dietary beta-carotene consumed (mcg per day).}
\item{\code{reedit}}{Dietary retinol consumed (mcg per day).}
\item{\code{betaplasma}}{Plasma beta-carotene (ng/ml).}
\item{\code{retplasma}}{Plasma Retinol (ng/ml).}
}
}

\source{\url{http://staff.pubhealth.ku.dk/~tag/Teaching/share/data/Plasma.html}}
\references{

Nierenberg, D. W., Stukel, T. A., Baron, J. A., Dain, B. J., and Greenberg, E. R. (1989).
	Determinants of plasma levels of beta-carotene and retinol. \emph{American Journal of Epidemiology}, \bold{130}, 511-521.

Meyer, M. C., Hackstadt, A. J., and Hoeting, J. A. (2011).
	Bayesian estimation and inference for generalized partial linear models using shape-restricted splines.
	\emph{Journal of Nonparametric Statistics}, \bold{23}(4), 867-884.

}
\examples{data(plasma)}

\keyword{datasets}
