\name{QGmvparams}
\alias{QGmvparams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Quantitative Genetics parameters from GLMM estimates (multivariate analysis).
}
\description{
  Function to estimate the multivariate quantitative genetics parameters on the observed scale based on estimates on the latent scale. Contrary to the univariate function, this one cannot use the analytical closed forms and yields a list of paramaters instead of a data.frame.
}
\usage{
    QGmvparams(mu = NULL, vcv.G, vcv.P, models, predict = NULL, rel.acc = 0.001, 
    width = 10, n.obs = NULL, theta = NULL, verbose = TRUE, mask = NULL)
}

\arguments{
  \item{mu}{Vector of latent intercepts estimated from a GLMM (ignored if predict is not \code{NULL}). (numeric)}
  \item{vcv.G}{Genetic additive variance-covariance matrix (a.k.a. G-matrix). (numeric)}
  \item{vcv.P}{Total phenotypic variance-covariance matrix. Usually, the sum of all the estimated variance-covariance matrices. (numeric)}
  \item{models}{A vector containing the names of the model used or a list which elements contain the list of the functions needed (inverse-link, distribution variance and derivative of the inverse-link, as stated in the output of \code{QGlink.funcs()}). (character vector or list of lists of functions) Available model names are :
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
  \item{rel.acc}{Relative accuracy of the integral approximation. (numeric)}
  \item{width}{Parameter for the integral computation. The default value is 10, which should be sensible for most models. (numeric)}
  \item{predict}{Optional matrix of predicted values on the latent scale (each trait in each column). The latent predicted values must be computed while only accounting for the fixed effects (marginal to the random effects). (numeric)}
  \item{n.obs}{Number of "trials" for the "binomN" distribution. (numeric)}
  \item{theta}{Dispersion parameter for the Negative Binomial distribution. The parameter \code{theta} should be such as the variance of the distribution is \code{mean + mean^2 / theta}. (numeric)}
  \item{verbose}{Should the function be verbose? (boolean)}
  \item{mask}{Masking filter for removing predictions that don't exist in the population (e.g. female predictions for males for a sex - based bivariate model). Should the same dimensions as \code{predict} and values should be \code{FALSE} when the predictions should be filtered out.}
}
\details{
The function typically uses integral numerical approximation provided by the R2Cuba package to compute multivariate quantitative genetics parameters on the observed scale, from latent estimates yielded by a GLMM. It cannot use closed form solutions.

Only the most typical distribution/link function couples are implemented through the \code{models} argument. If you used an "exotic" GLMM, you can provide a list containg lists of functions corresponding to the model. The list of functions should be implemented as is the output of \code{QGlink.funcs()}, i.e. three elements: the inverse link functions named \code{inv.link}, the derivative of this function named \code{d.inv.link} and the distribution variance named \code{var.func} (see Example below).

Some distributions require extra-arguments. This is the case for "binomN", which require the number of trials N, passed with the argument \code{n.obs}. The distribution "negbin" requires a dispersion parameter \code{theta}, such as the variance of the distribution is \code{mean + mean^2 / theta} (mean/dispersion parametrisation). For now, the arguments \code{n.obs} and \code{theta} can be used for ONE distribution only.

If fixed effects (apart from the intercept) have been included in the GLMM, they can be included through the argument \code{predict} as a matrix of the marginal predicted values, i.e. predicted values excluding the random effects, for each trait (one trait per column of the matrix, see Example below).Note this can considerably slow down the algorithm, especially when no closed form is used.
 
}

\value{
The function yields a list containing the following values:
  \item{mean.obs}{Vector of phenotypic means on the observed scale.}
  \item{vcv.P.obs}{Phenotypic variance-covariance matrix on the observed scale.}
  \item{vcv.G.obs}{Additive genetic variance-covariance (a.k.a. G-matrix) on the observed scale.}
}

\author{
Pierre de Villemereuil & Michael B. Morrissey
}
\seealso{
\code{\link{QGparams}}, \code{\link{QGlink.funcs}}, \code{\link{QGmvmean}}, \code{\link{QGvcov}}, \code{\link{QGmvpsi}}
}
\examples{
## Example using a bivariate model (Binary trait/Gaussian trait)
# Parameters
mu <- c(0, 1)
G <- diag(c(0.5, 2))
P <- diag(c(1, 4))

# Note: no phenotypic, nor genetic correlations, hence should be equal to univariate case!

# Multivariate analysis
QGmvparams(mu = mu, vcv.G = G, vcv.P = P, models = c("binom1.probit", "Gaussian"))
QGparams(mu = 0, var.a = 0.5, var.p = 1, model = "binom1.probit") # Consistent results!
# Reminder: the results are the same here because we have no correlation between the two traits

# Defining the model "by hand" using the list
list.models = list(
  model1 = list(inv.link = function(x){pnorm(x)},
                d.inv.link = function(x){dnorm(x)},
                var.func = function(x){pnorm(x) * (1 - pnorm(x))}),
  model2 = list(inv.link = function(x){x},
                d.inv.link = function(x){1},
                var.func = function(x){0})
)
# Running the same analysis than above
QGmvparams(mu = mu, vcv.G = G, vcv.P = P, models = list.models) # Same results!

# Using predicted values
# Say we have 100  individuals
n <- 100
# Let's simulate predicted values
p <- matrix(c(runif(n), runif(n)), ncol = 2)
# Note that p has as many as columns as we have traits (i.e. two)
# Multivariate analysis with predicted values
\dontrun{QGmvparams(predict = p, vcv.G = G, vcv.P = P, models = c("binom1.probit", "Gaussian"))}
}
