\name{cadmium}
\alias{cadmium}
\docType{data}
\title{
Cadmium dose-response meta data
}
\description{
This dataset includes minimal information of NCC-2012 meta data.
}
\usage{data("cadmium")}
\format{
  A data frame with 190 observations on the following 5 variables.
  \describe{
    \item{\code{gender}}{a numeric vector with 1 : Female, 0 : Male, 0.5 : Unknown or both}
    \item{\code{ethnicity}}{a integer vector with 1 : Asian, 2 : Caucasian}
    \item{\code{Ucd_GM}}{a numeric vector of Geometric means of urinary cadmium}
    \item{\code{b2_GM}}{a numeric vector of Geometric means of Beta2-Microglobulin}
    \item{\code{isOld}}{a logical vector whether the observation is older than 50}
  }
}
\references{
  Lee, Minjea, Choi, Taeryon; Kim, Jeongseon; Woo, Hae Dong (2013) Bayesian Analysis of Dose-Effect Relationship of Cadmium for Benchmark Dose Evaluation. \emph{Korean Journal of Applied Statistics}, \bold{26}(3), 453–470.

}
\examples{
\dontrun{
	data(cadmium)
}
}
\keyword{datasets}
