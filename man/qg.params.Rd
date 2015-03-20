\name{qg.params}
\alias{qg.params}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
qg.params(mu, var.a, var.p, model, width = 35, predict = NULL, closed.form = TRUE, custom.model = NULL, n.obs = NULL, theta = NULL, verbose = TRUE)
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
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{width}{
%%     ~~Describe \code{width} here~~
}
  \item{predict}{
%%     ~~Describe \code{predict} here~~
}
  \item{closed.form}{
%%     ~~Describe \code{closed.form} here~~
}
  \item{custom.model}{
%%     ~~Describe \code{custom.model} here~~
}
  \item{n.obs}{
%%     ~~Describe \code{n.obs} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (mu, var.a, var.p, model, width = 35, predict = NULL, 
    closed.form = TRUE, custom.model = NULL, n.obs = NULL, theta = NULL, 
    verbose = TRUE) 
{
    if (model == "binom1.probit" & closed.form) {
        if (verbose) 
            print("Using the closed forms for a Binomial1-probit model.")
        qg.binom1.probit(mu = mu, var.a = var.a, var.p = var.p, 
            predict = predict)
    }
    else if (model == "binomN.probit" & closed.form) {
        if (is.null(n.obs)) {
            stop("binomN.probit model used, but no observation number (n.obs) defined.")
        }
        if (verbose) 
            print("Using the closed forms for a BinomialN-probit model.")
        qg.binomN.probit(mu = mu, var.a = var.a, var.p = var.p, 
            predict = predict, n.obs = n.obs, width = width)
    }
    else if (model == "Poisson.log" & closed.form) {
        if (verbose) 
            print("Using the closed forms for a Poisson-log model.")
        qg.Poisson.log(mu = mu, var.a = var.a, var.p = var.p, 
            predict = predict)
    }
    else if (model == "Poisson.sqrt" & closed.form) {
        if (verbose) 
            print("Using the closed forms for a Poisson-sqrt model.")
        qg.Poisson.sqrt(mu = mu, var.a = var.a, var.p = var.p, 
            predict = predict)
    }
    else if (model == "negbin.log" & closed.form) {
        if (is.null(theta)) {
            stop("negbin model used, but theta not defined.")
        }
        if (verbose) 
            print("Using the closed forms for a NegativeBinomial-log model.")
        qg.negbin.log(mu = mu, var.a = var.a, var.p = var.p, 
            predict = predict, theta = theta)
    }
    else if (model == "negbin.sqrt" & closed.form) {
        if (is.null(theta)) {
            stop("negbin model used, but theta not defined.")
        }
        if (verbose) 
            print("Using the closed forms for a NegativeBinomial-sqrt model.")
        qg.negbin.sqrt(mu = mu, var.a = var.a, var.p = var.p, 
            predict = predict, theta = theta)
    }
    else {
        if (is.null(custom.model)) {
            funcs = link.funcs(model, n.obs = n.obs, theta = theta)
        }
        else {
            funcs = custom.model
        }
        if (verbose) 
            print("Computing observed mean...")
        y_bar = mean.obs(mu = mu, var = var.p, link.inv = funcs$inv.link, 
            width = width, predict = predict)
        if (verbose) 
            print("Computing variances...")
        var_exp = var.exp(mu = mu, var = var.p, link.inv = funcs$inv.link, 
            obs.mean = y_bar, width = width, predict = predict)
        var_dist = var.dist(mu = mu, var = var.p, var.func = funcs$var.func, 
            width = width, predict = predict)
        var_obs = var_exp + var_dist
        if (verbose) 
            print("Computing Psi...")
        Psi = psi(mu = mu, var = var.p, d.link.inv = funcs$d.inv.link, 
            width = width, predict = predict)
        data.frame(mean.obs = y_bar, var.obs = var_obs, Psi = Psi, 
            var.a.obs = (Psi^2) * var.a, h2.obs = ((Psi^2) * 
                var.a)/var_obs)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
