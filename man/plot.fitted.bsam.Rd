\name{plot.fitted.bsam}
\alias{plot.fitted.bsam}
\title{ Plot a fitted.bsam object}
\description{
Plots the data and the fit for Bayesian spectral analysis models.
}
\usage{
\method{plot}{fitted.bsam}(x, type, ask, ggplot2, legend.position, \dots)
}
\arguments{
  \item{x}{ a \code{fitted.bsam} object}
  \item{type}{ the type of fitted plot. The default is on the scale of the response variable\code{type="response"}; the alternative \code{type="term"} is on the scale of the nonparametric predictor. Note that this affects only on glm type models. For example, binomial model with the default option gives the predicted probabilites.}
  \item{ask}{see. \code{\link{par}}}
  \item{ggplot2}{a logical variable. If \code{TRUE} the \code{ggplot2} package is used.}
  \item{legend.position}{the position of legends (\dQuote{none}, \dQuote{left}, \dQuote{right}, \dQuote{bottom}, \dQuote{top}). It is used when \code{ggplot2 = TRUE}.}
  \item{\dots}{ other options to pass to the plotting functions}
}

\value{Returns a plot.}

\seealso{
\code{\link{bsaq}}, \code{\link{bsaqdpm}}, \code{\link{bsar}}, \code{\link{bsardpm}}, \code{\link{fitted.bsam}}
}

\examples{
## See examples for bsaq, bsaqdpm, bsar, and bsardpm
}
