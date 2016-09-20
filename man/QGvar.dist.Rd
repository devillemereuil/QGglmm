\name{QGvar.dist}
\alias{QGvar.dist}

\title{Compute the distribution variance}
\description{
This function computes the variance emerging from the error distribution around the individual expected value. This variance, added to the variance of the individual expected values themselves (see \code{\link{QGvar.exp}}) yields the total observed phenotypic variance.
}
\usage{
QGvar.dist(mu=NULL, var, var.func, predict = NULL, width = 10)
}
\arguments{
  \item{mu}{Latent intercept estimated from a GLMM (ignored if predict is not \code{NULL}). (numeric of length 1)}
  \item{var}{Latent total phenotypic variance estimated from a GLMM. Usually, the sum of the estimated variances of the random effects, plus the "residual" variance. (numeric of length 1)}
  \item{var.func}{Function giving the variance of the distribution according to a given latent value. (function)}
  \item{predict}{Optional vector of predicted values on the latent scale (i.e. matrix product \bold{Xb}). The latent predicted values must be computed while only accounting for the fixed effects (marginal to the random effects). (numeric)}
  \item{width}{Parameter for the integral computation. The integral is evaluated from \code{mu}-\code{width*sqrt(var)} to \code{mu}+\code{width*sqrt(var)}. The default value is 10, which should be sensible for most models. (numeric)}
}
\details{
The distribution variance is the part of the observed variance emerging from the error distribution. It is calculated as an average error variance over all possible latent values. The distribution variance added to the variance of the expected values gives the total phenotypic variance on the observed scale.

The variance function (\code{var.func}) is a function giving the variance of the error distribution of the GLMM according to a given latent value.

Using a Poisson distribution with a logarithm link, this function is \code{exp(x)}, because the variance of a Poisson is its mean. Using a Negative Binomial distribution with a logarithm link, this function will now be \code{exp(x)+exp(2*x)/theta}. Note that the dispersion parameter \code{theta} is necessary for a Negative Binomial distribution.

The \code{var.func} function is yielded by \code{\link{QGlink.funcs}} according to a given distribution.link model (see Example below).

Contrary to \code{\link{QGparams}}, \code{QGvar.exp} never uses the closed form solutions, but always compute the integrals.
}
\value{This function yields the distribution variance. (numeric)}

\author{
Pierre de Villemereuil & Michael B. Morrissey
}
\seealso{
\code{\link{QGvar.exp}}, \code{\link{QGparams}}, \code{\link{QGpred}}, \code{\link{QGlink.funcs}}, \code{\link{QGmean}}, \code{\link{QGpsi}}
}
\examples{
##Example using Poisson.log model
mu=1
va=0.2
vp=0.5

#The variance function is simply the inverse-link function
#because the variance of a Poisson is its mean
varfunc<-function(x){exp(x)}

QGvar.dist(mu=mu,var=vp,var.func=varfunc)

#The QGlink.funcs gives a ready-to-use var.func
funcs<-QGlink.funcs(name="Poisson.log")

#Calculating the distribution variance
vdist=QGvar.dist(mu=mu,var=vp,var.func=funcs$var.func)

vdist           #Same value as above

#Calculating the variance of the expected values
vexp=QGvar.exp(mu=mu,var=vp,link.inv=funcs$inv.link)

#The phenotypic variance on the observed scale is then:
vexp+vdist

#This computation is automatically performed by QGparams
#but directly using the closed form solutions when available
QGparams(mu=mu,var.p=vp,var.a=va,model="Poisson.log")
#var.obs is equal to the sum above
}

