\name{qg.negbin.sqrt}
\alias{qg.negbin.sqrt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
qg.negbin.sqrt(mu, var.a, var.p, theta, predict = NULL)
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
  \item{theta}{
%%     ~~Describe \code{theta} here~~
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
function (mu, var.a, var.p, theta, predict = NULL) 
{
    if (is.null(predict)) 
        predict = 0
    lambda = mean((mu + predict)^2 + var.p)
    var_obs = mean((mu + predict)^4 + 6 * var.p * ((mu + predict)^2) + 
        3 * (var.p^2)) - lambda^2 + lambda + (mean((mu + predict)^4 + 
        6 * var.p * ((mu + predict)^2) + 3 * (var.p^2))/theta)
    Psi = mean(2 * (mu + predict))
    data.frame(mean.obs = lambda, var.obs = var_obs, Psi = Psi, 
        var.a.obs = (Psi^2) * var.a, h2.obs = ((Psi^2) * var.a)/var_obs)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
