\name{QGrn_generate_shape}
\alias{QGrn_generate_shape}

\title{Compute the gradient of the reaction norm shape function}
\description{
Function yielding the gradient (vector of derivative)
}
\usage{
QGrn_generate_shape(expr, pars)
}
\arguments{
  \item{expr}{An expression giving the curve of the reaction norm, where the environment is named "x".}
  \item{pars}{A character vector with the names of all the parameters as used in \code{expr}.}
}
\details{
The function takes an expression (as often expressed for non-linear modelling frameworks) and returns a function that can be directly used by the \code{shape} argument of other \code{QGrn_*} functions.
}
\value{
A function to be passed on to the \code{shape} argument of \code{QGrn_vplas()} or \code{QGrn_vgen()}.
}

\author{
Pierre de Villemereuil
}

\seealso{
\code{\link{QGrn_vplas}}, \code{\link{QGrn_vgen}}, \code{\link{QGrn_mean_by_env}}
}
\examples{
# An expression of a reaction norm shape
gausgompz <- expression(
    cmax * exp(
        - exp(rho * (x - xopt) - 6) -
            sigmagaus * (x - xopt)^2
    ))

# Computing the function yielding the gradient
gausgompz <- QGrn_generate_shape(expr = gausgompz,
                                 pars = c("cmax", "xopt", "rho", "sigmagaus"))
}

