\name{QGmvmean}
\alias{QGmvmean}

\title{Compute the multivariate phenotypic mean on the observed scale}
\description{
This function calculates the multivariate phenotypic mean on the observed scale from multivariate latent mean and variance-covariance matrix.
}
\usage{
QGmvmean(mu=NULL, vcov, link.inv, predict = NULL, rel.acc=0.01, width = 10)
}
\arguments{
  \item{mu}{Vector of latent intercepts estimated from a GLMM (ignored if predict is not \code{NULL}). (numeric)}
  \item{vcov}{Latent total phenotypic variance-covariance matrix estimated from a GLMM. Usually, the sum of all the estimated variance-covariance matrices. (numeric)}
  \item{link.inv}{Inverse functions of the link functions. This function should accept a vector and yield a vector of the same length, see Details and Example below. (function)}
  \item{predict}{Optional matrix of predicted values on the latent scale (each trait in each column). The latent predicted values must be computed while only accounting for the fixed effects (marginal to the random effects). (numeric)}
  \item{rel.acc}{Relative accuracy of the integral approximation. (numeric)}
  \item{width}{Parameter for the integral computation. The default value is 10, which should be sensible for most models. (numeric)}
}

\details{
This function needs the multivariate latent population mean (\code{mu}) or the marginal predicted values (\code{predict}) and the total latent variance-covariance matrix (\code{vcov}) to compute the observed phenotypic mean.

To do so, it also requires the inverse functions of the link functions (\code{link.inv}). For an analysis with d traits, the function given to the \code{link.inv} argument should use a vector of length d and yield a vector of length d (see Example below).
}
\value{This function yields the mutlivariate phenotypic mean on the observed scale. (numeric)}

\author{
Pierre de Villemereuil & Michael B. Morrissey
}
\seealso{
\code{\link{QGmean}}, \code{\link{QGmvparams}}, \code{\link{QGlink.funcs}}, \code{\link{QGvcov}}, \code{\link{QGmvpsi}}
}
\examples{
##Example using a bivariate model (Binary trait/Gaussian trait)
#Parameters
mu=c(0,1)
P=diag(c(1,4))

#Note: no phenotypic, nor genetic correlations, hence should be equal to univariate case!

#Setting up the link functions
inv.links=function(vec){c(pnorm(vec[1]),vec[2])} #probit link, identity link respectively

#Computing the multivariate mean on observed scale
QGmvmean(mu=mu, vcov=P, link.inv=inv.links)
QGmean(mu=0,var=1,link.inv=pnorm)     #Same result than trait 1!
QGmean(mu=1,var=4,link.inv=function(x){x})     #Same result than trait 2!
#Reminder: the results are the same here because we have no correlation between the two traits
}
