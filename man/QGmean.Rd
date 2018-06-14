\name{QGmean}
\alias{QGmean}

\title{Compute the phenotypic mean on the observed scale}
\description{
This function calculates the phenotypic mean on the observed scale from the latent mean and variance.
}
\usage{
QGmean(mu = NULL, var, link.inv, predict = NULL, width = 10)
}
\arguments{
  \item{mu}{Latent intercept estimated from a GLMM (ignored if predict is not \code{NULL}). (numeric of length 1)}
  \item{var}{Latent total variance estimated from a GLMM. Usually, the sum of the estimated variances of the random effects, plus the "residual" variance. (numeric of length 1)}
  \item{link.inv}{Inverse function of the link function. (function)}
  \item{predict}{Optional vector of predicted values on the latent scale (i.e. matrix product \bold{Xb}). The latent predicted values must be computed while only accounting for the fixed effects (marginal to the random effects). (numeric)}
  \item{width}{Parameter for the integral computation. The integral is evaluated from \code{mu} - \code{width * sqrt(var)} to \code{mu} + \code{width * sqrt(var)}. The default value is 10, which should be sensible for most models. (numeric)}
}

\details{
This function needs the latent population mean (\code{mu}) or the marginal predicted values (\code{predict}) and the total latent variance (i.e. total latent variance \code{var}) to compute the observed phenotypic mean. To do so, it also requires the inverse function of the link function.
  
For example, if the link function is the natural logarithm, the inverse-link function will be the exponential. The inverse-link functions for many models are yielded by the \code{\link{QGlink.funcs}} function.

Contrary to \code{\link{QGparams}}, \code{QGmean.obs} never uses the closed form solutions, but always compute the integrals.
}
\value{This function yields the phenotypic mean on the observed scale. (numeric)}

\author{
Pierre de Villemereuil & Michael B. Morrissey
}
\seealso{
\code{\link{QGmvmean}}, \code{\link{QGparams}}, \code{\link{QGpred}}, \code{\link{QGlink.funcs}}, \code{\link{QGvar.dist}}, \code{\link{QGvar.exp}}, \code{\link{QGpsi}}
}
\examples{
## Computing the observed mean for a probit link
QGmean(mu = 0.3, var = 1, link.inv = pnorm)
# The theoretical expectation is
1 - pnorm(0, 0.3, sqrt(1 + 1))

# Or, using the QGlink.funcs function
QGmean(mu = 0.3, var = 1, link.inv = QGlink.funcs(name = "binom1.probit")$inv.link)

## Computing the observed mean for a logarithm link
QGmean(mu = 1, var = 1, link.inv = exp)
# The theoretical expectation is
exp(1 + 0.5 * 1)

# This computation is automatically performed by QGparams
# but directly using the closed form solution when available
QGparams(mu = 1, var.p = 1, var.a = 0.5, model = "Poisson.log")

}
