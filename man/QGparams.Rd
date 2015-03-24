\name{QGparams}
\alias{QGparams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Quantitative Genetics parameters from GLMM estimates.
}
\description{
  Function to estimate the quantitative genetics parameters on the observed scale based on estimates on the latent scale. The function yields a data.frame which includes the phenotypic mean and variance, as well as the additive genetic variance and heritability, on the observed scale.
}
\usage{
QGparams(mu, var.a, var.p, model="", width = 35, predict = NULL, closed.form = TRUE, custom.model = NULL, n.obs = NULL, theta = NULL, verbose = TRUE)
}

\arguments{
  \item{mu}{Latent intercept estimated from a GLMM. (numeric)}
  \item{var.a}{Latent additive genetic variance estimated from a GLMM. (numeric)}
  \item{var.p}{Latent total phenotypic variance estimated from a GLMM. Usually, the sum of the estimated variances of the random effects, plus the "residual" variance. (numeric)}
  \item{model}{Name of the used model, i.e. distribution.link. Ignored if \code{custom.model} is not NULL. (character) Available models are :
    \itemize{
      \item{"binom1.probit"} Binomial with 1 trial (binary data) with a probit link
      \item{"binomN.probit"} Binomial with N tria with a probit link (require the parameter \code{n.obs})
      \item{"binom1.logit"} Binomial with 1 trial (binary) with a logit link
      \item{"binomN.logit"} Binomial with N trial with a logit link (require the parameter \code{n.obs})
      \item{"Poisson.log"} Poisson distribution wiht a log link
      \item{"Poisson.sqrt"} Poisson distribution with a square-root link
      \item{"negbin.log"} Negative-Binomial distribution wiht a log link (require the parameter \code{theta})
      \item{"negbin.sqrt"} Negative-Binomial distribution with a square-root link (require the parameter \code{theta})
    }
  }
  \item{width}{Parameter for the integral computation. The integral is evaluated from \code{mu}-\code{width*sqrt(var.p)} to \code{mu}+\code{width*sqrt(var.p)}. The default value is 35, which should be sensible for most models. (numeric)}
  \item{predict}{Optional vector of predicted values on the latent scale (i.e. matrix product \bold{Xb}). The latent predicted values must be computed while only accounting for the fixed effects (marginal to the random effects). (numeric)}
  \item{closed.form}{When available, should closed forms be used instead of integral computations? (boolean)}
  \item{custom.model}{If the model used is not available using the \code{model} argument, a list of functions describing the model can be provided. (list of functions, see Details)}
  \item{n.obs}{Number of "trials" for the "binomN" distribution. (numeric)}
  \item{theta}{Dispersion parameter for the Negative Binomial distribution. The parameter \code{theta} should be such as the variance of the distribution is \code{mean + mean**2/theta}. (numeric)}
  \item{verbose}{Should the function be verbose? (boolean)}
}
\details{
The function typically uses precise integral numerical approximation to compute quantitative genetics parameters on the observed scale, from latent estimates yielded by a GLMM. If closed form solutions for the integrals are available, it uses them if \code{closed.form=TRUE}.

Only the most typical distribution/link function couples are implemented in the function. If you used an "exotic" GLMM, you can use the \code{custom.model} argument. It should take the form of a list of functions. The first function should be the inverse of the link function named \code{inv.link}, the second function should be the "distribution variance" function named \code{var.func} and the third function should be the derivative of the inverse link function named \code{d.inv.link} (see Example below).

Some distributions require extra-arguments. This is the case for "binomN", which require the number of trials N, passed with the argument \code{n.obs}. The distribution "negbin" requires a dispersion parameter \code{theta}, such as the variance of the distribution is \code{mean + mean**2/theta} (mean/dispersion parametrisation).
}

\value{
The function yields a data.frame containing the following values:
  \item{mean.obs}{Phenotypic mean on the observed scale.}
  \item{var.obs}{Phenotypic variance on the observed scale.}
  \item{var.a.obs}{Additive genetic variance on the observed scale.}
  \item{h2.obs}{Heritability on the observed scale.}
}
\references{
\strong{REF!}
}

\author{
Pierre de Villemereuil & Micheal B. Morrissey
}

\examples{
##Example using binary data
#Parameters
mu=0
va=1
vp=2

#Simulating data l = mu + a + e
lat<-mu+rnorm(1000,0,sqrt(va))+rnorm(1000,0,sqrt(vp-va))
y<-rbinom(1000,1,pnorm(lat))

#Expected results
QGparams(mu=0,var.p=2,var.a=1,model="binom1.probit")
#Simulated results for mean and variance
mean(y)
var(y)

#Using integral approximations
QGparams(mu=0,var.p=2,var.a=1,model="binom1.probit",closed.form=FALSE)
#Note that the approximation is exactly equal to the results obtained with the closed form

#Let's create a custom model
custom=list(inv.link=function(x){pnorm(x)},
    var.func=function(x){pnorm(x)*(1-pnorm(x))},
    d.inv.link=function(x){dnorm(x)})
    
QGparams(mu=0,var.p=2,var.a=1,custom.model=custom)
}