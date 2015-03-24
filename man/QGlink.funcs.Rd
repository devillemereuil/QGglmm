\name{QGlink.funcs}
\alias{QGlink.funcs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
List of functions according to a distribution and a link function
}
\description{
Function yielding different functions (inverse-link, variance function, derivative of the inverse-link) according to a distribution and link function.
}
\usage{
link.funcs(name, n.obs = NULL, theta = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{Name of the distribution.link couple. See }
  \item{n.obs}{
%%     ~~Describe \code{n.obs} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
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
\code{\link{QGparams}}, \code{\link{QGpredict}}, \code{\link{QGmean.obs}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (link.name, n.obs = NULL, theta = NULL) 
{
    if (link.name == "binom1.probit") {
        inv.link = function(x) {
            pnorm(x)
        }
        var.func = function(x) {
            pnorm(x) * (1 - pnorm(x))
        }
        d.inv.link = function(x) {
            dnorm(x)
        }
    }
    else if (link.name == "binomN.probit") {
        if (is.null(n.obs)) {
            stop("binomN model used, but no observation number (n.obs) defined.")
        }
        inv.link = function(x) {
            n.obs * pnorm(x)
        }
        var.func = function(x) {
            n.obs * pnorm(x) * (1 - pnorm(x))
        }
        d.inv.link = function(x) {
            n.obs * dnorm(x)
        }
    }
    else if (link.name == "binom.logit") {
        inv.link = function(x) {
            plogis(x)
        }
        var.func = function(x) {
            plogis(x) * (1 - plogis(x))
        }
        d.inv.link = function(x) {
            dlogis(x)
        }
    }
    else if (link.name == "binomN.logit") {
        if (is.null(n.obs)) {
            stop("binomN model used, but no observation number (n.obs) defined.")
        }
        inv.link = function(x) {
            n.obs * plogis(x)
        }
        var.func = function(x) {
            n.obs * plogis(x) * (1 - plogis(x))
        }
        d.inv.link = function(x) {
            n.obs * dlogis(x)
        }
    }
    else if (link.name == "threshold") {
        stop("Not implemented yet")
    }
    else if (link.name == "Poisson.log") {
        inv.link = function(x) {
            exp(x)
        }
        var.func = function(x) {
            exp(x)
        }
        d.inv.link = function(x) {
            exp(x)
        }
    }
    else if (link.name == "Poisson.sqrt") {
        inv.link = function(x) {
            x^2
        }
        var.func = function(x) {
            x^2
        }
        d.inv.link = function(x) {
            2 * x
        }
    }
    else if (link.name == "negbin.log") {
        if (is.null(theta)) {
            stop("negbin model used, but theta not defined.")
        }
        inv.link = function(x) {
            exp(x)
        }
        var.func = function(x) {
            exp(x) + (exp(2 * x)/theta)
        }
        d.inv.link = function(x) {
            exp(x)
        }
    }
    else if (link.name == "negbin.sqrt") {
        if (is.null(theta)) {
            stop("negbin model used, but theta not defined.")
        }
        inv.link = function(x) {
            x^2
        }
        var.func = function(x) {
            (x^2) + ((x^4)/theta)
        }
        d.inv.link = function(x) {
            2 * x
        }
    }
    else {
        stop("Invalid model name. Use a valid model name or enter a custom model specification.")
    }
    list(inv.link = inv.link, var.func = var.func, d.inv.link = d.inv.link)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
