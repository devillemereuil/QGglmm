\name{QGicc}
\alias{QGicc}

\title{
  Intra - Class Correlation coefficients (ICC) on the observed data scale
}
\description{
  Function to estimate the Intra - Class Correlation coefficients (ICC, a.k.a. repeatability - like estimates) on the observed scale based on estimates on the latent scale. For a specific variance component, the function yields a data.frame which includes the phenotypic mean and variance, as well as the variance component and associated ICC, on the observed data scale.
}
\usage{
    QGicc(mu = NULL, var.comp, var.p, model = "", width = 10, predict = NULL,
    closed.form = TRUE, custom.model = NULL, n.obs = NULL, theta = NULL, verbose = TRUE)
}

\arguments{
  \item{mu}{Latent intercept estimated from a GLMM (ignored if predict is not \code{NULL}). (numeric of length 1)}
  \item{var.comp}{Latent variance component for which ICC needs to be computed, estimated from a GLMM. (numeric of length 1)}
  \item{var.p}{Latent total phenotypic variance estimated from a GLMM. Usually, the sum of the estimated variances of the random effects, plus the "residual" variance. (numeric of length 1)}
  \item{model}{Name of the used model, i.e. distribution.link. Ignored if \code{custom.model} is not NULL. (character) Available models are :
    \itemize{
      \item{"Gaussian"} Gaussian distribution with identity link (e.g. LMM)
      \item{"binom1.probit"} Binomial with 1 trial (binary data) with a probit link
      \item{"binomN.probit"} Binomial with N tria with a probit link (require the parameter \code{n.obs})
      \item{"binom1.logit"} Binomial with 1 trial (binary) with a logit link
      \item{"binomN.logit"} Binomial with N trial with a logit link (require the parameter \code{n.obs})
      \item{"Poisson.log"} Poisson distribution wiht a log link
      \item{"Poisson.sqrt"} Poisson distribution with a square - root link
      \item{"negbin.log"} Negative - Binomial distribution wiht a log link (require the parameter \code{theta})
      \item{"negbin.sqrt"} Negative - Binomial distribution with a square - root link (require the parameter \code{theta})
    }
  }
  \item{width}{Parameter for the integral computation. The integral is evaluated from -\code{width * sqrt(var.comp)} to \code{width * sqrt(var.comp)}. The default value is 10, which should be sensible for most models. (numeric)}
  \item{predict}{Optional vector of predicted values on the latent scale (i.e. matrix product \bold{Xb}). The latent predicted values must be computed while only accounting for the fixed effects (marginal to the random effects). (numeric)}
  \item{closed.form}{When available, should closed forms be used instead of integral computations? (boolean)}
  \item{custom.model}{If the model used is not available using the \code{model} argument, a list of functions describing the model can be provided. (list of functions, see Details)}
  \item{n.obs}{Number of "trials" for the "binomN" distribution. (numeric)}
   \item{theta}{Dispersion parameter for the Negative Binomial distribution. The parameter \code{theta} should be such as the variance of the distribution is \code{mean + mean^2 / theta}. (numeric)}
  \item{verbose}{Should the function be verbose? (boolean)}
}
\details{
The function typically uses precise integral numerical approximation to compute parameters on the observed scale, from latent estimates yielded by a GLMM. If closed form solutions for the integrals are available, it uses them if \code{closed.form = TRUE}.

Only the most typical distribution/link function couples are implemented in the function. If you used an "exotic" GLMM, you can use the \code{custom.model} argument. It should take the form of a list of functions. The first function should be the inverse of the link function named \code{inv.link}, the second function should be the "distribution variance" function named \code{var.func} and the third function should be the derivative of the inverse link function named \code{d.inv.link} (see Example below).

Some distributions require extra-arguments. This is the case for "binomN", which require the number of trials N, passed with the argument \code{n.obs}. The distribution "negbin" requires a dispersion parameter \code{theta}, such as the variance of the distribution is \code{mean + mean^2 / theta} (mean/dispersion parametrisation).

If fixed effects (apart from the intercept) have been included in the GLMM, they can be included as marginal predicted values, i.e. predicted values excluding the random effects, which can be calculated as the matrix product \strong{Xb} where \strong{X} is the design matrix and \strong{b} is the vector of fixed effects estimates. To do so, provide the vector of marginal predicted values using the argument \code{predict}. Note this can considerably slow down the algorithm, especially when no closed form is used.
}

\value{
The function yields a data.frame containing the following values:
  \item{mean.obs}{Phenotypic mean on the observed scale.}
  \item{var.obs}{Phenotypic variance on the observed scale.}
  \item{var.comp.obs}{Component variance on the observed scale.}
  \item{icc.obs}{ICC on the observed scale.}
}

\author{
Pierre de Villemereuil & Michael B. Morrissey
}
\seealso{
\code{\link{QGparams}}, \code{\link{QGpred}}, \code{\link{QGlink.funcs}}, \code{\link{QGmean}}, \code{\link{QGvar.dist}}, \code{\link{QGvar.exp}}, \code{\link{QGpsi}}
}
\examples{
## Example using Poisson count data
# Parameters
mu <- 0
va <- 0.5
vm <- 0.2  # Maternal effect
vp <- 1

# Simulating data l = mu + a + e
lat <- mu + 
       rnorm(1000, 0, sqrt(va)) + 
       rnorm(1000, 0, sqrt(vm)) +
       rnorm(1000, 0, sqrt(vp - (va + vm)))
y   <- rpois(1000, exp(lat))

# Computing the broad - sense heritability
QGicc(mu = mu, var.p = vp, var.comp = va, model = "Poisson.log")
# Computing the maternal effect ICC
QGicc(mu = mu, var.p = vp, var.comp = vm, model = "Poisson.log")

# Using integral computation
QGicc(mu = mu, var.p = vp, var.comp = vm, model = "Poisson.log", closed.form = FALSE)
# Note that the "approximation" is exactly equal to the results obtained with the closed form

# Let's create a custom model
custom <- list(inv.link = function(x){exp(x)},
    var.func = function(x){exp(x)},
    d.inv.link = function(x){exp(x)})
    
QGicc(mu = mu, var.p = vp, var.comp = vm, custom.model = custom)
# Again, exactly equal

# Integrating over a posterior distribution
# e.g. output from MCMCglmm named "model"
# df <- data.frame(mu = model$Sol[, 'intercept'], 
#                  vm = model$VCV[, 'mother'], 
#                  vp = rowSums(model$VCV))
# params <- apply(df, 1, function(row){
#        QGicc(mu = row$mu, var.comp = row$vm, var.p = row$vp, model = "Poisson.log")
# })
}
