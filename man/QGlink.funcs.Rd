\name{QGlink.funcs}
\alias{QGlink.funcs}

\title{List of functions according to a distribution and a link function}
\description{
Function yielding different functions (inverse-link, variance function, derivative of the inverse-link) according to a distribution and link function.
}
\usage{
QGlink.funcs(name, n.obs = NULL, theta = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{Name of the distribution.link couple. See \code{\link{QGparams}} for a complete list of model available. (character)}
  \item{n.obs}{Optional parameter required for "binomN" distributions (number of "trials"). See \code{\link{QGparams}}. (numeric)}
  \item{theta}{Optional parameter required for "negbin" distributions (dispersion parameter). See \code{\link{QGparams}}. (numeric)}
}
\details{
This function takes the name of a distribution.link couple and yields several important functions such as the inverse-link function and its derivative, as well as the "distribution variance function".

The inverse-link function is the inverse function of the link function. For example, if the link function is the natural logarithm (typically for a Poisson distribution), then the inverse-link function is the exponential.

The distribution variance function is a function yielding the variance of the distribution for a given latent trait. For a Poisson distribution, the variance is equal to the mean, hence the variance function is equal to the inverse-link function. For a binomial distribution, the variance is N*p(l)*(1-p(l)), where p is the inverse-link function.

For some distributions, such as "binomN" and "negbin", some extra-parameters are required.
}
\value{
This function yields a list of function:
  \itemize{
    \item{\code{inv.link}}{Inverse function of the link function. (function)}
    \item{\code{var.func}}{Distribution variance function. (function)}
    \item{\code{inv.link}}{Derivative of the inverse-link function. (function)}
  }
}

\author{
Pierre de Villemereuil & Micheal B. Morrissey
}

\seealso{
\code{\link{QGparams}}, \code{\link{QGmvparams}}
}
\examples{
##Getting the functions for a Poisson.log model
QGlink.funcs("Poisson.log")
#Note that because the variance is equal to the mean in a Poisson distribution
#and the derivative of exp is exp
#all functions are the same!

##Getting the functions for a binom1.probit model
QGlink.funcs("binom1.probit")

##The function QGparams automatically computes these functions
QGparams(mu=0,var.p=2,var.a=1,model="binom1.logit")
#Hence this is the same as using the custom.model argument with QGlink.funcs
QGparams(mu=0,var.p=2,var.a=1,custom.model=QGlink.funcs("binom1.logit"))

##We can create our own custom set of functions
#Let's create a custom model exactly identical to QGlink.funcs("binom1.logit")
custom=list(inv.link=function(x){plogis(x)},
    var.func=function(x){plogis(x)*(1-plogis(x))},
    d.inv.link=function(x){dlogis(x)})
    
QGparams(mu=0,var.p=2,var.a=1,custom.model=custom)
}

