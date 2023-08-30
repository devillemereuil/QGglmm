\name{QGrn_vplas}
\alias{QGrn_vplas}

\title{Compute the average plastic variance}
\description{
This function calculates the plastic variance of the average (over genotypes) reaction norm shape.
}
\usage{
QGrn_vplas(env, shape, theta, G_theta, width = 10, fixed = NA, correction = FALSE)
}
\arguments{
  \item{env}{Vector of environmental values (numeric).}
  \item{shape}{Function providing the shape of the reaction norm, the environment should correspond to the variable \code{x} and the parameters to a matrix \code{pars}, each row for each parameter (function).}
  \item{theta}{Average (in the genetic sense) parameters of the shape function. (numeric)}
  \item{G_theta}{Genetic variance-covariance matrix of the parameters. It can be of lesser dimensions than \code{theta}, see \code{fixed} parameter. (numerical matrix)}
  \item{width}{Parameter for the integral computation. The integral is evaluated from \code{mu} - \code{width * sqrt(var)} to \code{mu} + \code{width * sqrt(var)}. The default value is 10, which should be sensible for most models. (numeric)}
  \item{fixed}{If some parameters of \code{shape}, included in \code{theta} are not included in the \code{G_theta} matrix, then those dimensions are considered as genetically "fixed". Hence, \code{fixed} should contain a vector of the \bold{index} of those parameters. Otherwise (if all parameters vary genetically), and by default, fixed should be set as NA. (integer)}
  \item{correction}{Should Bessel's correction (dividing the sum-of-squares by N-1 rather than N) be used (TRUE) or not (FALSE, default). The default is FALSE, because it is likely that other components such as the total phenotypic variance is computed over the number of individuals (generally large number) rather than the number of environments (generally small number). The best is to manually use Bessel's correction over the proper number of data points. (boolean)}
} 

\details{
It is very important that the parameters are in the same order in the argument of the \code{shape} function, in \code{theta} and in \code{G_theta}.
}
\value{This function yields the plastic variance arising from the average shape of the reaction norm (numeric).}

\author{
Pierre de Villemereuil
}
\seealso{
\code{\link{QGrn_vgen}}, \code{\link{QGrn_va}}
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
gausgompz <- QGrn_generate_shape(expr, c("cmax", "xopt", "rho", "sigmagaus"))

# Theta
theta <- c(cmax = 1, xopt = 0.9, rho = 8, sigmagaus = 0.4)
# G, only for cmax and xopt
G     <- matrix(c(0.1,      0.01,
                  0.01,     0.05),
                ncol = 2)

# Computing V_plas
QGrn_vplas(env     = vec_env,
           shape   = gausgompz,
           theta   = theta,
           G_theta = G,
           width   = 10,
           fixed   = c(3, 4),
           correction = FALSE)
# Note that fixed is set for the third and forth parameters than are not in G

# Note the quite large difference when Bessel's correction is used
QGrn_vplas(env     = vec_env,
           shape   = gausgompz,
           theta   = theta,
           G_theta = G,
           width   = 10,
           fixed   = c(3, 4),
           correction = TRUE)
}
