\name{QGmvpsi}
\alias{QGmvpsi}

\title{Compute a multivariate "Psi" (used to compute the additive genetic variance on the observed scale).}
\description{
This function computes a multivariate version of the parameter "Psi" which relates the additive genetic variance-covariance matrix on the latent scale to the additive genetic variance-covariance matrix on the observed scale: G.obs = Psi \%*\% G \%*\% t(Psi)
}
\usage{
QGmvpsi(mu = NULL, vcov, d.link.inv, predict = NULL, 
        rel.acc = 0.001, width = 10, compound = NULL, mask = NULL)
}

\arguments{
  \item{mu}{Vector of latent intercepts estimated from a GLMM (ignored if predict is not \code{NULL}). (numeric)}
  \item{vcov}{Latent total phenotypic variance-covariance matrix estimated from a GLMM. Usually, the sum of all the estimated variance-covariance matrices. (numeric)}
  \item{d.link.inv}{Derivative of the inverse-link functions. This function should accept a vector and yield a vector of the same length, see Details and Example below. (function)}
  \item{predict}{Optional matrix of predicted values on the latent scale (each trait in each column). The latent predicted values must be computed while only accounting for the fixed effects (marginal to the random effects). (numeric)}
  \item{rel.acc}{Relative accuracy of the integral approximation. (numeric)}
  \item{width}{Parameter for the integral computation. The default value is 10, which should be sensible for most models. (numeric)}
  \item{compound}{A vector of two indices, or list of such vectors (e.g. \code{list(c(1,2), c(4,5))}), providing the locations of "compound" distributions in the input (i.e. the dimensions that need to be "merged" into one in the output). The input must be ordered so that the first component (1 and 4 in the example above) is the latent trait for positive values and the second (2 and 5 in the example above) is the latent trait for the zero-component. (integer of length 2, or list of such vectors)}
  \item{mask}{Masking filter for removing predictions that don't exist in the population (e.g. female predictions for males for a sex-based bivariate model). Should the same dimensions as \code{predict} and values should be \code{FALSE} when the predictions should be filtered out.}
}
\details{
The multivariate parameter "Psi" is a diagonal matrix which elements are the average of the derivative of the inverse-link function. The additive genetic variance-covariance matrix on the latent scale G is linked to the additive genetic variance-covariance matrix on the observed scale G.obs through Psi: G.obs = Psi \%*\% G \%*\% t(Psi).

This function requires the derivatives of the inverse-link functions (\code{d.link.inv}). For an analysis with d traits, the function given to the \code{d.link.inv} argument should use a vector of length d and yield a vector of length d (see Example below). When using \code{compound}, the function corresponding to the compound distribution should accept a 2-rows input and yield a 2-rows output.
}
\value{This function yields the matrix "Psi". (numeric)}

\author{
Pierre de Villemereuil & Michael B. Morrissey
}
\seealso{
\code{\link{QGpsi}}, \code{\link{QGmvparams}}, \code{\link{QGlink.funcs}}, \code{\link{QGvcov}}, \code{\link{QGmvpsi}}
}
\examples{
## Example using a bivariate model (Binary trait/Gaussian trait)
# Parameters
mu <- c(0, 1)
G <- diag(c(0.5, 2))
P <- diag(c(1, 4))

# Setting up the derivatives of the inverse-link functions
dinvs <- function(mat) {matrix(c(dnorm(mat[1, ]), rep(1, length(mat[2, ]))),
                               nrow = 2, 
                               byrow = TRUE)}
# The derivative of pnorm() is dnorm(), and the derivative of the identity is 1

# Computing Psi
Psi <- QGmvpsi(mu = mu, vcov = P, d.link.inv = dinvs)
# Computing genetic additive variance-covariance matrix on the observed scale
Psi
G.obs <- Psi \%*\% G \%*\% t(Psi)

QGparams(mu = 0, var.a = 0.5, var.p = 1, model = "binom1.probit")
# Same additive variance than trait 1
# Reminder: the results are the same here because we have no correlation between the two traits
}
