\name{fs}
\alias{fs}
\title{Specify a Fourier Basis Fit in a BSAM Formula}
\description{
A symbolic wrapper to indicate a nonparametric term in a formula argument to bsaq, bsaqdpm, bsar, bsardpm, and gbsar.
}
\usage{
fs(x)
}
\arguments{
  \item{x}{a vector of the univariate covariate for nonparametric component}
}


\examples{
\dontrun{

# fit x using a Fourier basis
y ~ w + fs(x)

# fit x1 and x2 using a Fourier basis
y ~ fs(x1) + fs(x2)

}
}