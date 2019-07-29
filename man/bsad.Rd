\name{bsad}
\alias{bsad}

\title{Bayesian Semiparametric Density Estimation}
\description{
	This function fits a semiparametric model, which consists of parametric and
	nonparametric components, for estimating density using a logistic Gaussian process.
}

\usage{

bsad(x, xmin, xmax, nint, MaxNCos, mcmc = list(), prior = list(),
smoother = c('geometric', 'algebraic'),
parametric = c('none', 'normal', 'gamma', 'laplace'), marginal.likelihood = TRUE,
verbose = FALSE)
}

\arguments{
  \item{x}{			a vector giving the data from which the density estimate is to be computed.}

  \item{xmin}{			minimum value of x. }

  \item{xmax}{			maximum value of x. }

  \item{nint}{			number of grid points for plots (need to be odd).
				The default is 201.}

  \item{MaxNCos}{		maximum number of Fourier coefficients. }

  \item{mcmc}{			a list giving the MCMC parameters.
				The list includes the following integers (with default values in parentheses):
                     		\code{kappaloop (5)} giving the number of MCMC loops within each choice of kappa,
				\code{nblow (10000)} giving the number of MCMC in transition period,
				\code{nskip (10)} giving the thinning interval,
				\code{smcmc (1000)} giving the number of MCMC for analysis, and
                     		\code{ndisp (1000)} giving the number of saved draws to be displayed on screen
				(the function reports on the screen when every
				\code{ndisp} iterations have been carried out). }

  \item{prior}{			a list giving the prior information. The list includes the following parameters
				(default values specify the non-informative prior):
				\code{gmax} giving maximum value for gamma (default = 5),
				\code{PriorProbs} giving prior probability of parametric and semiparametric models,
				\code{beta_m0} and \code{beta_v0} giving the hyperparameters
				for prior distribution of the parametric coefficients,
                     		\code{r0} and \code{s0} giving the hyperparameters of \eqn{\sigma^2} for the logits,
				\code{u0} and \code{v0} giving the hyperparameters of \eqn{\tau^2} for Fourier coefficients,
				\code{PriorKappa} and \code{KappaGrid} giving prior on the number of cosine terms.}

  \item{smoother}{		types of smoothing priors for Fourier coefficients. See Details.}

  \item{parametric}{		specifying a distribution of the parametric part to be test.}

  \item{marginal.likelihood}{	  a logical variable indicating whether the log marginal likelihood is calculated.}

  \item{verbose}{   a logical variable. If \code{TRUE}, the iteration number and the Metropolis acceptance rate are printed to the screen.}

}

\details{
  This generic function fits a semiparametric model, which consists of parametric and nonparametric, for density estimation (Lenk, 2003):
  \deqn{f(x | \beta, Z) = \frac{\exp[h(x)^\top\beta + Z(x)]}{\int_\mathcal{X} \exp[h(y)^\top\beta + Z(y)]dG(y)}}
  where \eqn{Z} is a zero mean, second-order Gaussian process with bounded, continuous covariance function. i.e.,
  \deqn{E[Z(x), Z(y)] = \sigma(x,y), \quad \int_\mathcal{X}ZdG = 0 ~~(a.s.)}

  Using the Karhunen-Loeve Expansion, \eqn{Z} is represented as infinite series with random coefficients
  \deqn{Z(x) = \sum_{j=1}^\infty \theta_j\varphi_j(x), }
  where \eqn{\{\varphi_j\}} is the cosine basis, \eqn{\varphi_j(x)=\sqrt{2}\cos[j\pi G(x)]}.

  For the random Fourier coefficients of the expansion, two smoother priors are assumed (optional),
  \deqn{\theta_j | \tau, \gamma \sim N(0, \tau^2\exp[-j\gamma]), ~ j \ge 1 ~ (geometric ~smoother)}
  \deqn{\theta_j | \tau, \gamma \sim N(0, \tau^2\exp[-ln(j+1)\gamma]), ~ j \ge 1 ~ (algebraic ~smoother)}

  The coefficient \eqn{\beta} have the popular normal prior,
  \deqn{\beta | m_{0,\beta}, V_{0,\beta} \sim N(m_{0,\beta}, V_{0,\beta})}

  To complete the model specification, independent hyper priors are assumed,
  \deqn{\tau^2 | r_0, s_0 \sim IGa(r_0/2, s_0/2)}
  \deqn{\gamma | w_0 \sim Exp(w_0)}

  Note that the posterior algorithm is based on computing a discrete version of the likelihood over a fine mesh on \eqn{\mathcal{X}}.
}

\value{
  An object of class \code{bsad} representing the Bayesian spectral analysis density estimation model fit.
  Generic functions such as \code{print}, \code{fitted} and \code{plot} have methods to show the results of the fit.

  The MCMC samples of the parameters in the model are stored in the list \code{mcmc.draws},
  the posterior samples of the fitted values are stored in the list \code{fit.draws}, and
  the MCMC samples for the log marginal likelihood are saved in the list \code{loglik.draws}.
  The output list also includes the following objects:

  \item{post.est}{	posterior estimates for all parameters in the model. }

  \item{lmarg}{		log marginal likelihood. }

  \item{ProbProbs}{	posterior probability of models. }

  \item{call}{		the matched call. }

  \item{mcmctime}{     running time of Markov chain from \code{system.time()}. }

}


\references{

Jo, S., Choi, T., Park, B. and Lenk, P. (2019). bsamGP: An R Package for Bayesian Spectral Analysis Models
 Using Gaussian Process Priors. \emph{Journal of Statistical Software}, \bold{90}, 310-320.


Lenk, P. (2003) Bayesian semiparametric density estimation and model verification
	using a logistic Gaussian process. \emph{Journal of Computational and Graphical Statistics}, \bold{12}, 548-565.

}

\examples{
\dontshow{
set.seed(1)
x <- rgamma(50,2,1)

# mcmc parameters
mcmc <- list(nblow = 100,
             smcmc = 100,
             nskip = 1,
             ndisp = 100,
             kappaloop = 1)

# fits BSAD model
fout <- bsad(x = x, xmin = 0, xmax = 10, nint = 501, mcmc = mcmc,
            smoother = 'geometric', parametric = 'gamma')

# Summary
print(fout); summary(fout)

# fitted values
fit <- fitted(fout, HPD = FALSE)

# predictive density plot
plot(fit, ask = FALSE)
}
\dontrun{
############################
# Old Faithful geyser data #
############################
data(faithful)
attach(faithful)

# mcmc parameters
mcmc <- list(nblow = 10000,
	           smcmc = 1000,
	           nskip = 10,
	           ndisp = 1000,
	           kappaloop = 5)

# fits BSAD model
fout <- bsad(x = eruptions, xmin = 0, xmax = 8, nint = 501, mcmc = mcmc,
             smoother = 'geometric', parametric = 'gamma')

# Summary
print(fout); summary(fout)

# fitted values
fit <- fitted(fout)

# predictive density plot
plot(fit, ask = TRUE)

detach(faithful)
}
}

\keyword{ logistic Gaussian process }
\keyword{ density estimation }
