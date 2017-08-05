\name{function.shape}
\alias{function.shape}

\title{
Shape of Nonlinear Function
}

\description{
This function returns types of functions with shape constraints.
}

\usage{
function.shape(shape = c("Free", "Increasing", "Decreasing", 
"IncreasingConvex", "DecreasingConcave", "IncreasingConcave", 
"DecreasingConvex", "IncreasingS", "DecreasingS", 
"IncreasingRotatedS", "DecreasingRotatedS", "InvertedU", "Ushape"))
}


\arguments{
  \item{shape}{ a vector giving types of shape restriction.}
}


\note{
  This function is called by \code{bsaq}, \code{bsaqdpm}, \code{bsar}, 
  \code{bsardpm}, and \code{gbsar}, not called by user.
}

\keyword{ shape-restricted regression }
