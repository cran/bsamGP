\name{London.Mortality}
\alias{London.Mortality}
\docType{data}
\title{Daily Moratlity in London}
\description{
  The London.Mortality data consists of daily death occurrences from Jan. 1st, 1993 to Dec. 31st, 2006 and corresponding weather observations including temperature and humidity in London.
}
\usage{data(London.Mortality)}
\format{
  A data frame with 5113 observations on the following 7 variables.
  \describe{
    \item{date}{date in YYYY-MM-DD.}
    \item{tmean}{Mean temperature.}
    \item{tmin}{Minimum dry-bulb temperature.}
    \item{tmax}{Maximum dry-bulb temperature.}
    \item{dewp}{Dew point.}
    \item{rh}{Relative humidity.}
    \item{death}{the number of death occurences.}
  }
}
\source{
	Office for National Statistics

	British Atmospheric Data Centre

 	https://github.com/gasparrini/2015_gasparrini_Lancet_Rcodedata
}
\references{
  Armstrong BG, Chalabi Z, Fenn B, Hajat S, Kovats S, Milojevic A, Wilkinson P (2011). Association of mortality with high temperatures in a temperate climate: England and Wales. \emph{Journal of Epidemiology & Community Health}, \bold{65}(4), 340–345.


  Gasparrini A, Armstrong B, Kovats S, Wilkinson P (2012). The effect of high temperatures on cause-specific mortality in England and Wales. \emph{Occupational and Environmental Medicine}, \bold{69}(1), 56–61.


  Gasparrini A, Guo Y, Hashizume M, Lavigne E, Zanobetti A, Schwartz J, Tobias A, Tong S, Rocklöv J, Forsberg B, et al.(2015). Mortality risk attributable to high and low ambient temperature: a multicountry observational study. 
	\emph{The Lancet}, \bold{386}(9991), 369-375. 

}
\examples{
\dontrun{
	data(London.Mortality)
}
}
\keyword{datasets}
