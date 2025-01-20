\name{QGlink.funcs}
\alias{QGlink.funcs}

\title{List of functions according to a distribution and a link function}
\description{
Function yielding different functions (inverse-link, variance function, derivative of the inverse-link) according to a distribution and link function.
}
\usage{
QGlink.funcs(name, n.obs = NULL, theta = NULL)
}
\arguments{
  \item{name}{Name of the distribution.link couple. (character) Available models are :
    \itemize{
      \item{"Gaussian"} Gaussian distribution with identity link (e.g. LMM)
      \item{"binom1.probit"} Binomial with 1 trial (binary data) with a probit link
      \item{"binomN.probit"} Binomial with N tria with a probit link (require the parameter \code{n.obs})
      \item{"binom1.logit"} Binomial with 1 trial (binary) with a logit link
      \item{"binomN.logit"} Binomial with N trial with a logit link (require the parameter \code{n.obs})
      \item{"Poisson.log"} Poisson distribution with a log link
      \item{"Poisson.sqrt"} Poisson distribution with a square-root link
      \item{"negbin.log"} Negative - Binomial distribution wiht a log link (require the parameter \code{theta})
      \item{"negbin.sqrt"} Negative - Binomial distribution with a square-root link (require the parameter \code{theta})
      \item{"ZIPoisson.log.logit"} Zero-inflated Poisson distribution with a log link for the Poisson part and logit link for the ZI part (compound distribution)
      \item{"HuPoisson.log.logit"} Hurdle Poisson distribution with a log link for the Poisson part and logit link for the hurdle part (compound distribution)
      \item{"ordinal"} Multiple threshold model for ordinal categorical traits (require the parameter \code{cut.points})
    }}
  \item{n.obs}{Optional parameter required for "binomN" distributions (number of "trials"). See \code{\link{QGparams}}. (numeric)}
  \item{theta}{Optional parameter required for "negbin" distributions (dispersion parameter). See \code{\link{QGparams}}. (numeric)}
}
\details{
This function takes the name of a distribution.link couple and yields several important functions such as the inverse-link function and its derivative, as well as the "distribution variance function".

The inverse-link function is the inverse function of the link function. For example, if the link function is the natural logarithm (typically for a Poisson distribution), then the inverse-link function is the exponential.

The distribution variance function is a function yielding the variance of the distribution for a given latent trait. For a Poisson distribution, the variance is equal to the mean, hence the variance function is equal to the inverse-link function. For a binomial distribution, the variance is N * p(l) * (1 - p(l)), where p is the inverse-link function.

For some distributions, such as "binomN" and "negbin", some extra-parameters are required. The "compound" distributions (for now, zero-inflated and hurdle Poisson) functions take a 2-rows input instead of a vector (because they have two latent traits).

}
\value{
This function yields a list of function:
  \item{inv.link}{Inverse function of the link function. (function)}
  \item{var.func}{Distribution variance function. (function)}
  \item{inv.link}{Derivative of the inverse-link function. (function)}
}

\author{
Pierre de Villemereuil & Michael B. Morrissey
}

\seealso{
\code{\link{QGparams}}, \code{\link{QGmvparams}}
}
\examples{
## Getting the functions for a Poisson.log model
QGlink.funcs("Poisson.log")
# Note that because the variance is equal to the mean in a Poisson distribution
# and the derivative of exp is exp
# all functions are the same!

## Getting the functions for a binom1.probit model
QGlink.funcs("binom1.probit")

## The function QGparams automatically computes these functions
QGparams(mu = 0, var.p = 2, var.a = 1, model = "binom1.logit")
# Hence this is the same as using the custom.model argument with QGlink.funcs
QGparams(mu = 0, var.p = 2, var.a = 1, custom.model = QGlink.funcs("binom1.logit"))

## We can create our own custom set of functions
# Let's create a custom model exactly identical to QGlink.funcs("binom1.logit")
custom <- list(inv.link = function(x){plogis(x)},
               var.func = function(x){plogis(x) * (1 - plogis(x))},
               d.inv.link = function(x){dlogis(x)})
    
QGparams(mu = 0, var.p = 2, var.a = 1, custom.model = custom)
}

