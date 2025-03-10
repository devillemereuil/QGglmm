\name{QGvcov}
\alias{QGvcov}

\title{Compute the phenotypic variance-covariance matrix on the observed / expected scale}
\description{
This function computes the total phenotypic variance-covariance matrix on the observed or expected scales.
}
\usage{
QGvcov(mu = NULL, vcov, link.inv, var.func, mvmean.obs = NULL,
       predict = NULL, rel.acc = 0.001, width = 10,
       exp.scale = FALSE, compound = NULL, mask = NULL)
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
  \item{compound}{A vector of two indices, or list of such vectors (e.g. \code{list(c(1,2), c(4,5))}), providing the locations of "compound" distributions in the input (i.e. the dimensions that need to be "merged" into one in the output). The input must be ordered so that the first component (1 and 4 in the example above) is the latent trait for positive values and the second (2 and 5 in the example above) is the latent trait for the zero-component. (integer of length 2, or list of such vectors)}
  \item{mask}{Masking filter for removing predictions that don't exist in the population (e.g. female predictions for males for a sex-based bivariate model). Should the same dimensions as \code{predict} and values should be \code{FALSE} when the predictions should be filtered out.}
}
\details{
This function needs the multivariate latent population mean (\code{mu}) or the marginal predicted values (\code{predict}) and the total latent variance-covariance matrix (\code{vcov}) to compute the phenotypic variance-covariance matrix on the observed scale (or on the expected scale if \code{exp.scale} is \code{TRUE}).

To do so, it also requires the inverse functions of the link functions (\code{link.inv}) and the distribution variance functions (\code{var.func}). For an analysis with d traits, the function given to these arguments should use a vector of length d and yield a vector of length d (see Example below). When using \code{compound}, the functions corresponding to the compound distribution(s) should accept a 2-rows input and yield a vector.
}
\value{This function yields the phenotypic variance-covariance on the observed or expected scale. (numeric)}

\author{
Pierre de Villemereuil & Michael B. Morrissey
}
\seealso{
\code{\link{QGvar.exp}}, \code{\link{QGvar.dist}}, \code{\link{QGmvparams}}, \code{\link{QGlink.funcs}}, \code{\link{QGmvpsi}}
}
\examples{
## Example using a bivariate model (Binary trait/Gaussian trait)
# Parameters
mu <- c(0, 1)
P <- diag(c(1, 4))

# Note: no phenotypic, nor genetic correlations, hence should be equal to univariate case!

# Setting up the link functions
# Note that since the use of "cubature" to compute the integrals,
# the functions must use a matrix as input and yield a matrix as output,
# each row corresponding to a trait
inv.links <- function(mat) {matrix(c(pnorm(mat[1, ]), mat[2, ]), nrow = 2, byrow = TRUE)}

# Setting up the distribution variance functions
var.funcs <- function(mat) {matrix(c(pnorm(mat[1, ]) * (1 - pnorm(mat[1, ])), 0 * mat[2, ]),
                                   nrow = 2, 
                                   byrow = TRUE)}
# The first row is p * (1 - p) (variance of a binomial)
# The second row is 0 because no extra distribution is assumed for a Gaussian trait

# Computing the multivariate mean on observed scale
# Phenotypic VCV matrix on observed scale
QGvcov(mu = mu, vcov = P, link.inv = inv.links, var.func = var.funcs) 
# Phenotypic VCV matrix on the expected scale
QGvcov(mu = mu, vcov = P, link.inv = inv.links, var.func = var.funcs, exp.scale = TRUE) 

QGvar.exp(mu = 0, var = 1, link.inv = pnorm) # Same variance on the expected scale
QGvar.exp(mu = 0, var = 1, link.inv = pnorm) +
    QGvar.dist(mu = 0, var = 1, var.func = function(x){pnorm(x) * (1 - pnorm(x))})
# Same variance on the observed scale
# Reminder: the results are the same here because we have no correlation between the two traits
}

