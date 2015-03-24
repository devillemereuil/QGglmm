\name{QGpred}
\alias{QGpred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict the evolutionary response to selection on the observed scale
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
QGpred(mu, var.a, var.p, model, fitness.func, width = 35, predict = NULL, custom.model = NULL, n.obs = NULL, theta = NULL, verbose = TRUE)
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
  \item{fitness.func}{
%%     ~~Describe \code{fitness.func} here~~
}
  \item{width}{
%%     ~~Describe \code{width} here~~
}
  \item{predict}{
%%     ~~Describe \code{predict} here~~
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
function (mu, var.a, var.p, model, fitness.func, width = 35, 
    predict = NULL, custom.model = NULL, n.obs = NULL, theta = NULL, 
    verbose = TRUE) 
{
    if (is.null(predict)) 
        predict = 0
    print("Computing mean fitness...")
    Wbar = mean(sapply(predict, function(pred_i) {
        integrate(f = function(x) {
            lat.fitness(x) * dnorm(x, mu + pred_i, sqrt(var.p))
        }, lower = mu + pred_i - width * sqrt(var.p), upper = mu + 
            pred_i + width * sqrt(var.p))$value
    }))
    print("Computing the latent selection...")
    sel = (mean(sapply(predict, function(pred_i) {
        integrate(f = function(x) {
            x * lat.fitness(x) * dnorm(x, mu + pred_i, sqrt(var.p))
        }, lower = mu + pred_i - width * sqrt(var.p), upper = mu + 
            pred_i + width * sqrt(var.p))$value
    })) - (mu * Wbar))/Wbar
    print("Computing the latent response... (this might take a while if predict is large)")
    exp_W_a <- function(vec) {
        sapply(vec, function(a) {
            mean(sapply(predict, function(pred_i) {
                integrate(f = function(x) {
                  lat.fitness(x) * dnorm(x, mu + a + pred_i, 
                    sqrt(var.p - var.a))
                }, lower = mu + a + pred_i - width * sqrt(var.p - 
                  var.a), upper = mu + a + pred_i + width * sqrt(var.p - 
                  var.a))$value
            }))
        })
    }
    resp = (integrate(f = function(a) {
        a * exp_W_a(a) * dnorm(a, 0, sqrt(var.a))
    }, lower = -width * sqrt(var.a), upper = width * sqrt(var.a))$value)/Wbar
    data.frame(mean.lat.fit = Wbar, pred.lat.sel = sel, pred.lat.resp = resp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
