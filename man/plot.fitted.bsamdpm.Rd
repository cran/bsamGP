\name{plot.fitted.bsamdpm}
\alias{plot.fitted.bsamdpm}
\title{ Plot a fitted.bsamdpm object}
\description{
Plots the data and the fit for Bayesian spectral analysis models with Dirichlet process mixture error.
}
\usage{
\method{plot}{fitted.bsamdpm}(x, ask, ggplot2, legend.position, \dots)
}
\arguments{
  \item{x}{ a \code{fitted.bsamdpm} object}
  \item{ask}{see. \code{\link{par}}}
  \item{ggplot2}{a logical variable. If \code{TRUE} the \code{ggplot2} package is used.}
  \item{legend.position}{the position of legends (\dQuote{none}, \dQuote{left}, \dQuote{right}, \dQuote{bottom}, \dQuote{top}). It is used when \code{ggplot2 = TRUE}.}
  \item{\dots}{ other options to pass to the plotting functions}
}

\value{Returns a plot.}

\seealso{
\code{\link{bsaqdpm}}, \code{\link{bsardpm}}, \code{\link{fitted.bsamdpm}}
}

\examples{
## See examples for bsaqdpm and bsardpm
}
