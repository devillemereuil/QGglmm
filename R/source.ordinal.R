#----------------------------------------------------------------
#	QGglmm R package (functions for ordinal variables)
#	Functions to estimate QG parameters from GLMM estimates
#----------------------------------------------------------------
#	Pierre de Villemereuil (2015)
#----------------------------------------------------------------

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more detailsls.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

##--------------------------------------Ordinal model

qg.ordinal <- function(mu = NULL, 
                       var.a,
                       var.p,
                       cut.points,
                       predict = NULL)
{  
    # Various checks
    if (length(mu) > 1 | length(var.a) != 1 | length(var.p) != 1) {
        stop("The parameters mu, var.a and var.p must be of length 1, 
             please check your input.")
    }
    
    if (is.null(predict)) {predict <- mu};
    
    if (cut.points[1] != -Inf | cut.points[length(cut.points)] != Inf) {
        stop("Please ensure to input all cut points, including -Inf and Inf")
    }
    
    # Get number of categories from cut-points
    nb.cat <- length(cut.points) - 1
    
    # Observed means
    p <- numeric(nb.cat)
    for (i in 1:nb.cat) {
        p[i] <- mean(pnorm(cut.points[i+1], predict, sqrt(var.p+1)) - 
                     pnorm(cut.points[i], predict, sqrt(var.p+1)))
    }
    
    # Observed variance-covariance
    Sigma <- matrix(0, nrow = nb.cat, ncol = nb.cat)
    for (i in 1:nb.cat) {
        for (j in 1:nb.cat) {
            if (i==j) {
                Sigma[i,i] <- p[i] * (1 - p[i])
            } else  {
                Sigma[i,j] <- - p[i] * p[j]
            }
        }
    }
    
    # Psi
    Psi <- numeric(nb.cat)
    for (i in 1:nb.cat) {
        Psi[i] <- mean(dnorm(cut.points[i], predict, sqrt(var.p+1)) - 
                       dnorm(cut.points[i+1], predict, sqrt(var.p+1)))
    }
    
    # Additive genetic variance-covariance matrix
    G <- Psi %*% t(Psi) * var.a
    
    # Heritabilities on the observed data scale
    h2 <- diag(G) / diag(Sigma)
    
    # Return a list of QG parameters on the observed scale
    list(mean.obs = p, vcv.P.obs = Sigma, vcv.G.obs = G, h2.obs = h2)
}
