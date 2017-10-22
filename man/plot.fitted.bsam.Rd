\name{plot.fitted.bsam}
\alias{plot.fitted.bsam}
\title{ Plot a fitted.bsam object}
\description{
Plots the data and the fit for Bayesian spectral analysis models.
}
\usage{
\method{plot}{fitted.bsam}(x, ask, ggplot2, \dots)
}
\arguments{
  \item{x}{ a \code{fitted.bsam} object}
  \item{ask}{see. \code{\link{par}}}
  \item{ggplot2}{a logical variable. If \code{TRUE} the \code{ggplot2} package is used.}
  \item{\dots}{ other options to pass to the plotting functions}
}

\value{Returns a plot.}

\seealso{
\code{\link{bsaq}}, \code{\link{bsaqdpm}}, \code{\link{bsar}}, \code{\link{bsardpm}}, \code{\link{fitted.bsam}}
}

\examples{
## See examples for bsaq, bsaqdpm, bsar, and bsardpm
}
