\name{QGvar.exp}
\alias{QGvar.exp}

\title{Compute the variance of expected values (i.e. the latent values after inverse-link transformation.)}
\description{This function computes the variance of the expected values, i.e. the variance of the latent values after transformation through the inverse-link function. This variance, added to the distribution variance, yields to the phenotypic variance on the observed scale.}
\usage{
QGvar.exp(mu=NULL, var, link.inv, obs.mean = NULL, predict = NULL, width = 10)
}

\arguments{
  \item{mu}{Latent intercept estimated from a GLMM (ignored if predict is not \code{NULL}). (numeric of length 1)}
  \item{var}{Latent total phenotypic variance estimated from a GLMM. Usually, the sum of the estimated variances of the random effects, plus the "residual" variance. (numeric of length 1)}
  \item{link.inv}{Inverse function of the link function. (function)}
  \item{obs.mean}{Optional parameter giving the phenotypic mean on the observed scale. Automatically computed if not provided. (numeric)}
  \item{predict}{Optional vector of predicted values on the latent scale (i.e. matrix product \bold{Xb}). The latent predicted values must be computed while only accounting for the fixed effects (marginal to the random effects). (numeric)}
  \item{width}{Parameter for the integral computation. The integral is evaluated from \code{mu}-\code{width*sqrt(var)} to \code{mu}+\code{width*sqrt(var)}. The default value is 10, which should be sensible for most models. (numeric)}
}
\details{
The variance of the expected values is the variance that directly arise from the variance of the latent values, but after transformation through the inverse-link function. For example, using a logarithm link, this is the variance of exp(l) where l is the latent trait.

To compute the variance, the function needs the phenotypic mean on the observed scale. If this was already computed, it can be provided using the optional argument \code{obs.mean}, which will save computing time. Otherwise (default), the function will compute the mean on the observed scale before computing the variance.

This variance, when added to the distribution variance (see \code{\link{QGvar.dist}}) yields the phenotypic variance on the observed scale.

The function required for \code{link.inv} is yielded by \code{\link{QGlink.funcs}} according to a given distribution.link model (see Example below).

Contrary to \code{\link{QGparams}}, \code{QGvar.dist} never uses the closed form solutions, but always compute the integrals.
}
\value{This function yields the variance of the expected values. (numeric)}

\author{
Pierre de Villemereuil & Michael B. Morrissey
}
\seealso{
\code{\link{QGvar.dist}}, \code{\link{QGparams}}, \code{\link{QGpred}}, \code{\link{QGlink.funcs}}, \code{\link{QGmean}}, \code{\link{QGpsi}}
}
\examples{
##Example using Poisson.log model
mu=1
va=0.2
vp=0.5

#The inverse-link for a logarithm link is the exponential
inv.link<-function(x){exp(x)}

#We can then calculate the variance of expected values
QGvar.exp(mu=mu,var=vp,link.inv=inv.link)

#The mean on the observed scale can be computed beforehand
y_bar=QGmean(mu=mu,var=vp,link.inv=inv.link)
QGvar.exp(mu=mu,var=vp,obs.mean=y_bar,link.inv=inv.link)

#The QGlink.funcs gives a ready-to-use inverse-link function
funcs<-QGlink.funcs(name="Poisson.log")

#Calculating the distribution variance
vexp=QGvar.exp(mu=mu,var=vp,obs.mean=y_bar,link.inv=funcs$var.func)

vexp           #Same value as above

#Calculating the associated distribution variance
vdist=QGvar.dist(mu=mu,var=vp,var.func=funcs$var.func)

#The phenotypic variance on the observed scale is then:
vexp+vdist

#This computation is automatically performed by QGparams
#but directly using the closed form solutions when available
QGparams(mu=mu,var.p=vp,var.a=va,model="Poisson.log")
#var.obs is equal to the sum above
}
