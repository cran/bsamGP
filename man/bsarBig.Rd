\name{bsarBig}
\alias{bsarBig}

\title{Bayesian Spectral Analysis Regression for Big data}
\description{
	This function fits a Bayesian spectral analysis regression model for Big data.
}
\usage{
bsarBig(formula, nbasis, nint, mcmc = list(), prior = list(), verbose = FALSE)
}

\arguments{
  \item{formula}{ an object of class \dQuote{\code{\link{formula}}} }
  \item{nbasis}{  number of cosine basis functions. }
  \item{nint}{    number of grid points where the unknown function is evaluated for plotting.
				The default is 500. }
  \item{mcmc}{    a list giving the MCMC parameters.
				The list includes the following integers (with default values in parentheses):
				\code{nblow (10000)} giving the number of MCMC in transition period,
				\code{nskip (10)} giving the thinning interval,
				\code{smcmc (1000)} giving the number of MCMC for analysis, and
        \code{ndisp (1000)} giving the number of saved draws to be displayed on screen
				(the function reports on the screen when every
				\code{ndisp} iterations have been carried out). }
  \item{prior}{   a list giving the prior information. The list includes the following parameters
				(default values specify the non-informative prior):
				\code{sigma2_m0} and \code{sigma2_v0} giving the prior mean and variance
				of the inverse gamma prior for the scale parameter of response,
				\code{tau2_m0}, \code{tau2_s0} and \code{w0} giving the prior mean and standard deviation of
				smoothing prior. }
  \item{verbose}{ a logical variable. If \code{TRUE}, the iteration number and the Metropolis acceptance rate are printed to the screen.}
}

\value{
  The MCMC samples of the parameters in the model are stored in the list \code{mcmc.draws} and
  the posterior samples of the fitted values are stored in the list \code{fit.draws}.
  The output list also includes the following objects:

  \item{post.est}{	posterior estimates for all parameters in the model. }

  \item{call}{		the matched call. }

  \item{mcmctime}{      running time of Markov chain from \code{system.time()}. }
}

\seealso{
\code{\link{bsar}}
}
\examples{
# Ttrue function
ftrue <- function(x){
  ft <- 7*exp(-3*x) + 2*exp(-70*(x-.6)^2) - 2 + 5*x
  return(ft)
}

# Generate data
set.seed(1)

nobs <- 100000 # Number of observations
sigmat <- .5 # True sigma
nxgrid <- 500 # number of grid points: approximate likelihood & plots

xdata <- runif(nobs) # Generate x values
fobst <- ftrue(xdata) # True f at observations
ydata <- fobst + sigmat*rnorm(nobs)

# Compute grid on 0 to 1
xdelta <- 1/nxgrid
xgrid <- seq(xdelta/2, 1-xdelta/2, xdelta)
xgrid <- matrix(xgrid,nxgrid)
fxgridt <- ftrue(xgrid) # True f on xgrid

# Fit data
fout <- bsarBig(ydata ~ xdata, nbasis = 50, nint = nxgrid, verbose = TRUE)

# Plots
smcmc <- fout$mcmc$smcmc
t <- 1:smcmc
par(mfrow=c(2,2))
matplot(t, fout$mcmc.draws$theta, type = "l", main = "Theta", xlab = "Iteration", ylab = "Draw")
plot(t, fout$mcmc.draws$sigma, type = "l", main = "Sigma", xlab = "Iteration", ylab = "Draw")
matplot(t, fout$mcmc.draws$tau, type = "l", main = "Tau", xlab = "Iteration", ylab = "Draw")
matplot(t, fout$mcmc.draws$gamma, type = "l", main = "Gamma", xlab = "Iteration", ylab = "Draw")

dev.new()
matplot(fout$fit.draws$xgrid, cbind(fxgridt, fout$post.est$fhatm, fout$post.est$fhatq),
        type = "l", main = "Regression Function", xlab = "X", ylab = "Y")

# Compute RMISE for regression function
sse <- (fout$post.est$fhatm - fxgridt)^2
rmise <- intgrat(sse, 1/nxgrid)
rmise <- sqrt(rmise)
rmise
}

\keyword{ Big Data }
\keyword{ Gaussian processes }
