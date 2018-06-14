\name{QGpred}
\alias{QGpred}

\title{Predict the evolutionary response to selection on the observed scale}
\description{
This function uses an assumed or measured fitness function to compute evolutionary response to selection on the observed scale. To do so a latent fitness function must be provided to the function. This fitness function is used to compute the evolutionary response on the latent scale.
}

\usage{
    QGpred(mu = NULL, var.a, var.p, fit.func, d.fit.func, width = 10, 
           predict = NULL, verbose = TRUE)
}

\arguments{
  \item{mu}{Latent intercept estimated from a GLMM (set to 0 if \code{predict} is not \code{NULL}). (numeric of length 1)}
  \item{var.a}{Latent additive genetic variance estimated from a GLMM. (numeric of length 1)}
  \item{var.p}{Latent total phenotypic variance estimated from a GLMM. Usually, the sum of the estimated variances of the random effects, plus the "residual" variance. (numeric of length 1)}
  \item{fit.func}{Function giving the expected fitness on the observed scale for a given latent trait (see Example). (function)}
  \item{d.fit.func}{Derivative of the expected fitness to the latent trait (see Example). (function)}
  \item{width}{Parameter for the integral computation. The integral is evaluated from \code{mu} - \code{width * sqrt(var.p)} to \code{mu} + \code{width * sqrt(var.p)}. The default value is 10, which should be sensible for most models. (numeric)}
  \item{predict}{Optional vector of predicted values on the latent scale (i.e. matrix product \bold{Xb}). The latent predicted values must be computed while only accounting for the fixed effects (marginal to the random effects). (numeric)}
  \item{verbose}{Should the function be verbose? (boolean)}
}
\details{
The function uses the latent fitness function (\code{fit.func}) and latent quantitative genetics parameters to compute the expected selection differential and response on the latent scale.

There is no argument to describe the model used as it is already and implicitely contained in the calculation of \code{fit.func}.

If fixed effects were included during the estimation of the quantitative genetics parameters,  they can be included as marginal predicted values, i.e. predicted values excluding the random effects, which can be calculated as the matrix product \strong{Xb} where \strong{X} is the design matrix and \strong{b} is the vector of fixed effects estimates. To do so, provide the vector of marginal predicted values using the argument \code{predict}. Note this will considerably slow down the algorithm.

The predictions can be transposed on the observed scale by using the \code{\link{QGmean}} function (see Example below).
}
\value{
The function yields a data.frame containing:
  \itemize{
    \item{\code{mean.lat.fitness}}{Average latent fitness. (numeric)}
    \item{\code{lat.grad}}{Latent selection gradient. (numeric)}
    \item{\code{lat.sel}}{Latent selection differential. (numeric)}
    \item{\code{lat.resp}}{Latent evolutionary response to selection. (numeric)}
  }
}

\author{
Pierre de Villemereuil & Michael B. Morrissey
}

\seealso{
\code{\link{QGparams}}, \code{\link{QGlink.funcs}}, \code{\link{QGmean}}, \code{\link{QGvar.dist}}, \code{\link{QGvar.exp}}, \code{\link{QGpsi}}
}

\examples{
## Example with binary traits and a fitness measurement
# Let's assume we dispose of a binary trait measurement 
# and associated fitness of trait 0 (say 1) and trait 1 (say 1.86)
# We further assume a GLMM with Binomial distribution and probit link with:
mu <- -0.1
va <- 2
vp <- 2.5  # note that the latent heritability is very high

# Creating the latent fitness function
# i.e. expected fitness given a latent trait l
# We have a trait 1 with probability pnorm(l) with fitness 1.86
# We have a trait 0 with probability (1 - pnorm(l)) with fitness 1
lat.fit<- function(l){(1 - pnorm(l)) * 1 + pnorm(l) * 1.86}
# Derivate of the fitnes function
d.lat.fit<- function(l){- dnorm(l) * 1 + dnorm(l) * 1.86}

# Predicting the latent evolutionary response
pred <- QGpred(mu = mu, var.p = vp, var.a = va, fit.func = lat.fit, d.fit.func = d.lat.fit)

# Predicting the observed evolutionary response
# Current observed phenotypic mean
QGmean(mu = mu, var = vp, link.inv = QGlink.funcs("binom1.probit")$inv.link)
# Predicted observed phenotypic mean after selection
QGmean(mu = mu + pred$lat.resp, var = vp, link.inv = QGlink.funcs("binom1.probit")$inv.link)
}

