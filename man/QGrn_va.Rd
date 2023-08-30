\name{QGrn_va}
\alias{QGrn_va}

\title{Compute the genetic variance arising from the reaction norm}
\description{
This function calculates the genetic variance arising from genetic variation in the reaction norm.
}
\usage{
QGrn_va(env, d_shape, theta, G_theta, width = 10, fixed = NA, average = TRUE)
}
\arguments{
  \item{env}{Vector of environmental values (numeric).}
  \item{d_shape}{Function providing the \bold{gradient} of shape of the reaction norm, the environment should correspond to the variable \code{x} and the parameters to a matrix \code{pars}, each row for each parameter. The function should output a matrix of a number of rows equal to the number of parameters against which the gradient is computed. (function).}
  \item{theta}{Average (in the genetic sense) parameters of the shape function. (numeric)}
  \item{G_theta}{Genetic variance-covariance matrix of the parameters. It can be of lesser dimensions than \code{theta}, see \code{fixed} parameter. (numerical matrix)}
  \item{width}{Parameter for the integral computation. The integral is evaluated from \code{mu} - \code{width * sqrt(var)} to \code{mu} + \code{width * sqrt(var)}. The default value is 10, which should be sensible for most models. (numeric)}
  \item{fixed}{If some parameters of \code{shape}, included in \code{theta} are not included in the \code{G_theta} matrix, then those dimensions are considered as genetically "fixed". Hence, \code{fixed} should contain a vector of the \bold{index} of those parameters. Otherwise (if all parameters vary genetically), and by default, fixed should be set as NA. (integer)}
  \item{average}{Should the values for each environment be given, or just their average (default)? (boolean)}
} 

\details{
It is very important that the parameters are in the same order in the argument of the \code{shape} function, in \code{theta} and in \code{G_theta}.
}
\value{
This function yields the additive genetic variance (V_A), along with the gamma-decomposition of the impact of each parameter (and their covariance). If average is FALSE, those values are provided for each environment and a column is added to provide the environmental values. (data.frame)
}

\author{
Pierre de Villemereuil
}
\seealso{
\code{\link{QGrn_vgen}}
}
\examples{
# Some environment vector
vec_env <- seq(-2, 2)

# Shape function
expr <- expression(
     cmax * exp(
         - exp(rho * (x - xopt) - 6) -
             sigmagaus * (x - xopt)^2
     ))
d_gausgompz <- QGrn_generate_gradient(expr, c("cmax", "xopt"), c("cmax", "xopt", "rho", "sigmagaus"))

# Theta
theta <- c(cmax = 1, xopt = 0.9, rho = 8, sigmagaus = 0.4)
# G, only for cmax and xopt
G     <- matrix(c(0.1,      0.01,
                  0.01,     0.05),
                ncol = 2)

# Computing V_A and the gamma-decomposition
QGrn_va(env     = vec_env,
        d_shape = d_gausgompz,
        theta   = theta,
        G_theta = G,
        width   = 10,
        fixed   = c(3, 4),
        average = TRUE)
# Note that fixed is set for the third and forth parameters than are not in G

# Computing V_A and the gamma-decomposition for each environment
QGrn_va(env     = vec_env,
        d_shape = d_gausgompz,
        theta   = theta,
        G_theta = G,
        width   = 10,
        fixed   = c(3, 4),
        average = FALSE)
}
