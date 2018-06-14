\name{QGpsi}
\alias{QGpsi}

\title{Compute "Psi" (used to compute the additive genetic variance on the observed scale).}
\description{
This function computes the parameter "Psi" which relates the additive genetic variance on the latent scale to the additive genetic variance on the observed scale: Va.obs = (Psi^2) * Va
}
\usage{
QGpsi(mu = NULL, var, d.link.inv, predict = NULL, width = 10)
}

\arguments{
  \item{mu}{Latent intercept estimated from a GLMM (set to 0 if \code{predict} is not \code{NULL}). (numeric of length 1)}
  \item{var}{Latent total phenotypic variance estimated from a GLMM. Usually, the sum of the estimated variances of the random effects, plus the "residual" variance. (numeric of length 1)}
  \item{d.link.inv}{Derivative of the inverse-link function. (function)}
  \item{predict}{Optional vector of predicted values on the latent scale (i.e. matrix product \bold{Xb}). The latent predicted values must be computed while only accounting for the fixed effects (marginal to the random effects). (numeric)}
  \item{width}{Parameter for the integral computation. The integral is evaluated from \code{mu} - \code{width * sqrt(var)} to \code{mu} + \code{width * sqrt(var)}. The default value is 10, which should be sensible for most models. (numeric)}
}
\details{
The parameter "Psi" is the average of the derivative of the inverse-link function. The additive genetic variance on the observed scale is linked to the additive genetic variance on the latent scale by : Va.obs = (Psi^2) * Va.lat.
}
\value{This function yields the "Psi" parameter. (numeric)}

\author{
Pierre de Villemereuil & Michael B. Morrissey
}
\seealso{
\code{\link{QGmvpsi}}, \code{\link{QGparams}}, \code{\link{QGpred}}, \code{\link{QGlink.funcs}}, \code{\link{QGmean}}, \code{\link{QGvar.dist}}, \code{\link{QGvar.exp}}
}
\examples{
## Example using binom1.probit model
mu <- 0
va <- 1
vp <- 2
# The inverse-link for a probit is the CDF of a standard Gaussian
# Hence its derivative is the PDF of a standard Gaussian
dinv <- function(x){dnorm(x)}

# Computing Psi
Psi <- QGpsi(mu = 0, var = 2, d.link.inv = dinv)
# Computing additive variance on the observed scale
(Psi^2) * va

# This function is used by QGparams to obtain var.a.obs
QGparams(mu = 0, var.p = vp, var.a = va, model = "binom1.probit")
# Same results as above!
}
