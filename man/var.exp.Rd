\name{var.exp}
\alias{var.exp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
var.exp(mu, var, link.inv, obs.mean = NULL, width = 35, predict = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{var}{
%%     ~~Describe \code{var} here~~
}
  \item{link.inv}{
%%     ~~Describe \code{link.inv} here~~
}
  \item{obs.mean}{
%%     ~~Describe \code{obs.mean} here~~
}
  \item{width}{
%%     ~~Describe \code{width} here~~
}
  \item{predict}{
%%     ~~Describe \code{predict} here~~
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
function (mu, var, link.inv, obs.mean = NULL, width = 35, predict = NULL) 
{
    if (is.null(obs.mean)) {
        obs.mean = mean.obs(mu = mu, var = var, link.inv = link.inv, 
            width = width, predict = predict)
    }
    if (is.null(predict)) 
        predict = 0
    mean(sapply(predict, function(pred_i) {
        integrate(f = function(x) {
            (link.inv(x)^2) * dnorm(x, mu + pred_i, sqrt(var))
        }, lower = mu + pred_i - width * sqrt(var), upper = mu + 
            pred_i + width * sqrt(var))$value
    })) - (obs.mean^2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
