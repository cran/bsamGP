\name{plot.fitted.bsad}
\alias{plot.fitted.bsad}
\title{Plot a fitted.bsad object}
\description{
Plots the predictive density for Bayesian density estimation model using logistic Gaussian process
}
\usage{
\method{plot}{fitted.bsad}(x, ggplot2, legend.position, nbins, \dots)
}
\arguments{
  \item{x}{ a \code{fitted.bsad} object}
  \item{ggplot2}{a logical variable. If \code{TRUE} the \code{ggplot2} package is used.}
  \item{legend.position}{the position of legends (\dQuote{none}, \dQuote{left}, \dQuote{right}, \dQuote{bottom}, \dQuote{top}). It is used when \code{ggplot2 = TRUE}.}
  \item{nbins}{Number of bins used. Default is 30.}
  \item{\dots}{ other options to pass to the plotting functions}
}

\value{Returns a plot.}

\seealso{
\code{\link{bsad}}, \code{\link{fitted.bsad}}
}

\examples{
## See example for bsad
}
