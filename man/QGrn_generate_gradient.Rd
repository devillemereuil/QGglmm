\name{QGrn_generate_gradient}
\alias{QGrn_generate_gradient}

\title{Compute the gradient of the reaction norm shape function}
\description{
Function yielding the gradient (vector of derivative)
}
\usage{
QGrn_generate_gradient(expr, dpars, allpars)
}
\arguments{
  \item{expr}{An expression giving the curve of the reaction norm, where the environment is named "x".}
  \item{dpars}{A character vector with the names of the vector against which the gradient must be calculated (must be the same as used in \code{expr}.}
  \item{allpars}{A character vector with the names of all the parameters as used in \code{expr}.}
}
\details{
The function internally uses the utilities for computer algebra in R, e.g. the \code{D()} function, to compute the gradient of the reaction norm curve. It might fail to find the gradient for complex functions, e.g. based on more than \code{+}, \code{-}, \code{*}, \code{/}, \code{exp()} and \code{log()}.
}
\value{
A function to be passed on to the \code{d_shape} argument of \code{QGrn_va()}.
}

\author{
Pierre de Villemereuil
}

\seealso{
\code{\link{QGrn_va}}
}
\examples{
# An expression of a reaction norm shape
gausgompz <- expression(
    cmax * exp(
        - exp(rho * (x - xopt) - 6) -
            sigmagaus * (x - xopt)^2
    ))

# Computing the function yielding the gradient
d_gausgompz <- QGrn_generate_gradient(expr = gausgompz,
                                        dpars = c("cmax", "xopt"),
                                        allpars = c("cmax", "xopt", "rho", "sigmagaus"))
}

