\name{QGvcov}
\alias{QGvcov}

\title{Compute the phenotypic variance-covariance matrix on the observed/expected scale}
\description{
This function computes the total phenotypic variance-covariance matrix on the observed or expected scales.
}
\usage{
QGvcov(mu, vcov, link.inv, var.func, mvmean.obs = NULL,
          predict = NULL, rel.acc = 0.01, width = 10, exp.scale = FALSE)
}
\arguments{
  \item{mu}{Vector of latent intercepts estimated from a GLMM (ignored if predict is not \code{NULL}). (numeric)}
  \item{vcov}{Latent total phenotypic variance-covariance matrix estimated from a GLMM. Usually, the sum of all the estimated variance-covariance matrices. (numeric)}
  \item{link.inv}{Inverse functions of the link functions. This function should accept a vector and yield a vector of the same length, see Details and Example below. (function)}
  \item{var.func}{Function giving the variance function for each trait. This function should accept a vector and yield a vector of the same length, see Details and Example below. (function)}
    \item{mvmean.obs}{Optional parameter giving the multivariate phenotypic mean on the observed scale. Automatically computed if not provided. (numeric)}
  \item{predict}{Optional matrix of predicted values on the latent scale (each trait in each column). The latent predicted values must be computed while only accounting for the fixed effects (marginal to the random effects). (numeric)}
  \item{rel.acc}{Relative accuracy of the integral approximation. (numeric)}
  \item{width}{Parameter for the integral computation. The default value is 10, which should be sensible for most models. (numeric)}
  \item{exp.scale}{Should the variance-covariance matrix be computed on the expected scale? \code{FALSE} by default, which means the variance-covariance matrix is computed on the observed scale. (boolan)}
}
\details{
This function needs the multivariate latent population mean (\code{mu}) or the marginal predicted values (\code{predict}) and the total latent variance-covariance matrix (\code{vcov}) to compute the phenotypic variance-covariance matrix on the observed scale (or on the expected scale if \code{exp.scale} is \code{TRUE}).

To do so, it also requires the inverse functions of the link functions (\code{link.inv}) and the distribution variance functions (\code{var.func}). For an analysis with d traits, the function given to these arguments should use a vector of length d and yield a vector of length d (see Example below).
}
\value{This function yields the phenotypic variance-covariance on the observed or expected scale. (numeric)}

\author{
Pierre de Villemereuil & Micheal B. Morrissey
}
\seealso{
\code{\link{QGvar.exp}}, \code{\link{QGvar.dist}}, \code{\link{QGmvparams}}, \code{\link{QGlink.funcs}}, \code{\link{QGmvpsi}}
}
\examples{
##Example using a bivariate model (Binary trait/Gaussian trait)
#Parameters
mu=c(0,1)
P=diag(c(1,4))

#Note: no phenotypic, nor genetic correlations, hence should be equal to univariate case!

#Setting up the link functions
inv.links=function(vec){c(pnorm(vec[1]),vec[2])} #probit link, identity link respectively

#Setting up the distribution variance functions
var.funcs=function(vec){c(pnorm(vec[1])*(1-pnorm(vec[1])),0)} #probit and identity link respectively
#The first element of the vector is p*(1-p) (variance of a binomial)
#The second element is 0 because the identity link does not create variance

#Computing the multivariate mean on observed scale
#Phenotypic VCV matrix on observed scale
QGvcov(mu=mu, vcov=P, link.inv=inv.links, var.func=var.funcs) 
#Phenotypic VCV matrix on the expected scale
QGvcov(mu=mu, vcov=P, link.inv=inv.links, var.func=var.funcs,exp.scale=TRUE) 

QGvar.exp(mu=0,var=1,link.inv=pnorm) #Same variance on the expected scale
QGvar.exp(mu=0,var=1,link.inv=pnorm)+
QGvar.dist(mu=0,var=1,var.func=function(x){pnorm(x)*(1-pnorm(x))})
#Same variance on the observed scale
#Reminder: the results are the same here because we have no correlation between the two traits
}

