\name{qg.binomN.probit}
\alias{qg.binomN.probit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
qg.binomN.probit(mu, var.a, var.p, n.obs, predict = NULL, width = 35)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{var.a}{
%%     ~~Describe \code{var.a} here~~
}
  \item{var.p}{
%%     ~~Describe \code{var.p} here~~
}
  \item{n.obs}{
%%     ~~Describe \code{n.obs} here~~
}
  \item{predict}{
%%     ~~Describe \code{predict} here~~
}
  \item{width}{
%%     ~~Describe \code{width} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (mu, var.a, var.p, n.obs, predict = NULL, width = 35) 
{
    if (is.null(predict)) 
        predict = 0
    p = n.obs * mean(1 - pnorm(0, mu + predict, sqrt(var.p + 
        1)))
    prob.sq.int = mean(sapply(predict, function(pred_i) {
        integrate(f = function(x) {
            (pnorm(x)^2) * dnorm(x, mu + pred_i, sqrt(var.p))
        }, lower = mu + pred_i - width * sqrt(var.p), upper = mu + 
            pred_i + width * sqrt(var.p))$value
    }))
    var_obs = ((n.obs^2) - n.obs) * prob.sq.int - p^2 + p
    Psi = n.obs * mean(dnorm(0, (mu + predict), sqrt(var.p + 
        1)))
    data.frame(mean.obs = p, var.obs = var_obs, Psi = Psi, var.a.obs = (Psi^2) * 
        var.a, h2.obs = ((Psi^2) * var.a)/var_obs)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
