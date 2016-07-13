\name{QGmvpred}
\alias{QGmvpred}

\title{Predict the evolutionary response to selection on the observed scale}
\description{
This function uses an assumed or measured fitness function to compute evolutionary response to selection on the observed scale. To do so a latent fitness function must be provided to the function. This fitness function is used to compute the evolutionary response on the latent scale.
}

\usage{
    QGmvpred(mu, vcv.G, vcv.P, fit.func, d.fit.func, predict=NULL, rel.acc=0.01, 
    width = 10, verbose = TRUE)
}

\arguments{
  \item{mu}{Vector of latent intercepts estimated from a GLMM (ignored if predict is not \code{NULL}). (numeric)}
  \item{vcv.G}{Genetic additive variance-covariance matrix (a.k.a. G-matrix). (numeric)}
  \item{vcv.P}{Total phenotypic variance-covariance matrix. Usually, the sum of all the estimated variance-covariance matrices. (numeric)}
  \item{fit.func}{Function giving the expected fitness on the observed scale for a given latent trait (see Example). (function)}
  \item{d.fit.func}{Derivative of the expected fitness to the latent trait. This function should return a vector containing the partial derivative to each trait (see Example). (function)}
  \item{rel.acc}{Relative accuracy of the integral approximation. (numeric)}
  \item{width}{Parameter for the integral computation. The default value is 10, which should be sensible for most models. (numeric)}
  \item{predict}{Optional matrix of predicted values on the latent scale (each trait in each column). The latent predicted values must be computed while only accounting for the fixed effects (marginal to the random effects). (numeric)}
  \item{verbose}{Should the function be verbose? (boolean)}
}
\details{
The function uses the latent fitness function (\code{fitness.func}) and latent quantitative genetics parameters to compute the expected selection differential and response on the latent scale.

There is no argument to describe the model used as it is already and implicitely contained in the calculation of \code{fit.func} and \code{d.fit.func} (see Example below).

If fixed effects were included during the estimation of the quantitative genetics parameters,  they can be included as marginal predicted values, i.e. predicted values excluding the random effects, which can be calculated as the matrix product \strong{Xb} where \strong{X} is the design matrix and \strong{b} is the vector of fixed effects estimates. To do so, provide the vector of marginal predicted values using the argument \code{predict}. Note this will considerably slow down the algorithm.

The predictions can be transposed on the observed scale by using the \code{\link{QGmvmean}} function (see Example below).
}
\value{
The function yields a data.frame containing:
  \itemize{
    \item{\code{mean.lat.fitness}}{Average latent fitness. (numeric)}
    \item{\code{lat.grad}}{Latent selection gradient. (numeric)}
    \item{\code{lat.sel}}{Latent selection differential. (numeric)}
    \item{\code{lat.resp}}{Latent evolutionary response to selection. (numeric)}
  }
}

\author{
Pierre de Villemereuil & Micheal B. Morrissey
}

\seealso{
\code{\link{QGparams}}, \code{\link{QGlink.funcs}}, \code{\link{QGmean}}, \code{\link{QGvar.dist}}, \code{\link{QGvar.exp}}, \code{\link{QGpsi}}
}

\examples{
##Bivariate example with a binary trait and a Gaussian one

#Assume a bivariate GLMM with Binomial(probit)/Gaussian distributions with:
mu=c(0,10)
G=matrix(c(0.5,0,0,1),nrow=2)
P=matrix(c(1,0,0,2),nrow=2) 

#Link functions
inv.links=function(vec){c(pnorm(vec[1]),vec[2])}

#Creating the expected fitness function
#i.e. expected fitness given a latent trait vector l
#Say if the binary trait is 1, then the fitness is 0.5*"the Gaussian trait"
#But if the binary trait is 0, then the fitness is 0
lat.fit<-function(l){pnorm(l[1])*0.5*l[2]}
#Derivative of the above function
#This function yields a vector which elements are the derivative according to each trait
d.lat.fit<-function(l){c(dnorm(l[1])*0.5*l[2],pnorm(l[1])*0.5)}

#Predicting the latent evolutionary response
pred<-QGmvpred(mu=mu,vcv.P=P,vcv.G=G,fit.func=lat.fit,d.fit.func=d.lat.fit)



#Predicting the observed evolutionary response
#Current observed phenotypic mean
QGmvmean(mu=mu,vcov=P,link.inv=inv.links)
#Predicted observed phenotypic mean after selection
QGmvmean(mu=mu+pred$lat.resp,vcov=P,link.inv=inv.links)
}

