\name{QGglmm-package}
\alias{QGglmm-package}
\alias{QGglmm}
\docType{package}
\title{
Quantitative Genetics parameters for Generalised Linear Mixed Models
}
\description{
This package gives the values on the observed scale for several quantitative genetics parameter using estimates from a Generalised Linear Mixed Model (GLMM). If a fitness function is assumed or measured, it also predicts the evolutionary response to selection on the observed scale.

The two main functions of this package are \code{\link{QGparams}} and \code{\link{QGpred}}. The first allows to compute the quantitative genetics parameters on the observed scale for any given GLMM and its estimates. The second allows to compute a predicted response to evolution on the observed scale using GLMM estimates and an assumed/measured/inferred fitness function.

For some distribution/link models (e.g. Binomial/probit and Poisson and Negative Binomial with logartihm or square-root link), a closed form solutions of the integrals computed by this package are available. They are automatially used by \code{\link{QGparams}} and this function only.
}
\details{
  \tabular{ll}{
  Package: \tab QGglmm\cr
  Type: \tab Package\cr
  Version: \tab 0.4\cr
  Date: \tab 2016-07-13\cr
  License: \tab GPL-2\cr
  }
}