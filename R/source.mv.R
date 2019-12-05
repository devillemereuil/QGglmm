#----------------------------------------------------------------
#	QGglmm R package (multivariate)
#	Functions to estimate QG parameters from GLMM estimates
#----------------------------------------------------------------
#	Pierre de Villemereuil & Michael Morrissey (2015)
#----------------------------------------------------------------

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -----------------------------Helper functions-------------------------- ##

# Log-determinant of a VCV matrix
calc_logdet <- function(Sigma) {
    sum(log(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values))
}

# Vectorised multivariate Gaussian density function
# Shamelessly stolen from cubature vignette (credit to Balasubramanian Narasimhan)
# logdet = sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
vec_mvnorm <- function(x, mean, Sigma, logdet = NULL) {
    # If logdet not provided, compute it
    if (is.null(logdet)) {
        logdet <- calc_logdet(Sigma)
    }
    # Compute Mahalanobis distance (corresponds to the exp. part of the density)
    distval <- stats::mahalanobis(t(x), center = mean, cov = Sigma)
    # Compute the vectorised MVN density
    out <- exp(matrix(-(nrow(x) * log(2 * pi) + logdet + distval)/2, ncol = ncol(x)))
    return(out)
}

# Vectorised function to get upper-triangle squared matrix
vec_sq_uptri <- function(x) {
    apply(x, 2, function(col) {
        mat <- (col) %*% t(col)
        return(mat[upper.tri(mat, diag = TRUE)])
    })
}

##  ----------------------------General functions------------------------- ##

#Calculating the observed/expected scale mean (multivariate)
QGmvmean <- function(mu = NULL,
                     vcov, 
                     link.inv, 
                     predict = NULL, 
                     rel.acc = 0.001, 
                     width = 10, 
                     mask = NULL) 
{
    #Setting the integral width according to vcov (lower mean-w, upper mean+w)
    w <- sqrt(diag(vcov)) * width
    
    #Number of dimensions
    d <- length(w)
    
    #If predict is not included, then use mu, and 
    if(is.null(predict)) {
        if(is.null(mu)) {
            stop("Please provide either mu or predict.")
        } else {
            predict <- matrix(mu, nrow = 1)
        }
    }
    
    # Computing the logdet of vcov
    logdet <- calc_logdet(vcov)
    
    #Computing the mean
    #The double apply is needed to compute the mean for all "predict" values,
    #then average over them
    mat <- apply(predict, 1, function(pred_i) {
        cubature::hcubature(
            f  = function(x) {
                link.inv(x) * matrix(rep(vec_mvnorm(x, pred_i, vcov, logdet), d),
                                     nrow = d,
                                     byrow = TRUE)
            },
            lowerLimit = pred_i - w,
            upperLimit = pred_i + w,
            fDim       = d,
            tol        = rel.acc,
            absError   = 0.0001,
            vectorInterface = TRUE
        )$integral
    })
    
    #Applyign the mask if provided
    if (!is.null(mask)) {mat[t(mask)] <- NA}
    
    return(apply(mat, 1, mean, na.rm = TRUE))
}

#Calculating the expected scale variance-covariance matrix
QGvcov <- function(mu = NULL, 
                   vcov,
                   link.inv,
                   var.func,
                   mvmean.obs = NULL,
                   predict = NULL,
                   rel.acc = 0.001, 
                   width = 10,
                   exp.scale = FALSE,
                   mask = NULL) 
{
    #Setting the integral width according to vcov (lower mean-w, upper mean+w)
    w <- sqrt(diag(vcov)) * width
    
    #Number of dimensions
    d <- length(w)
    
    #If no fixed effects were included in the model
    if(is.null(predict)) { 
        if(is.null(mu)) {
            stop("Please provide either mu or predict.")
        } else {
            predict <- matrix(mu, nrow = 1)
        }
    }
    
    # Computing the logdet of vcov
    logdet <- calc_logdet(vcov)
    
    #Computing the upper-triangle matrix of "expectancy of the square"
    v <- apply(predict, 1,
               function(pred_i) {
                   cubature::hcubature(
                       f  = function(x) {
                           vec_sq_uptri(link.inv(x)) * 
                               matrix(rep(vec_mvnorm(x, pred_i, vcov, logdet), (d^2 + d) / 2),
                                      nrow = (d^2 + d) / 2,
                                      byrow = TRUE)
                       },
                       lowerLimit = pred_i - w,
                       upperLimit = pred_i + w,
                       fDim       = (d^2 + d) / 2,
                       tol        = rel.acc,
                       absError   = 0.0001,
                       vectorInterface = TRUE
                   )$integral
               })
    
    #Applying the mask if provided
    if(!is.null(mask)) {
        mask2 <- matrix(FALSE, ncol = ncol(v), nrow = nrow(v))
        mask2[((1:d) * ((1:d) + 1)) / 2, ] <- t(mask)
        v[mask2] <- NA
    }
    v <- apply(v, 1, mean, na.rm = TRUE)
    
    #Creating the VCV matrix
    vcv <- matrix(NA, d, d)
    vcv[upper.tri(vcv, diag = TRUE)] <- v
    vcv[lower.tri(vcv)] <- vcv[upper.tri(vcv)]
    
    #If necessary, computing the observed multivariate mean
    if (is.null(mvmean.obs)) {
        mvmean.obs <- 
            QGmvmean(mu,
                     vcov,
                     link.inv,
                     predict = predict,
                     rel.acc = rel.acc,
                     width = width)
    }
    
    #Computing the VCV matrix using Keonig's formuka
    vcv <- vcv - mvmean.obs %*% t(mvmean.obs)
    
    #Adding the distribution variance if needed (if exp.scale == FALSE)
    if (!exp.scale) {
        vardist <- apply(predict, 1, function(pred_i) {
            cubature::hcubature(
                f  = function(x) {
                    var.func(x) * matrix(rep(vec_mvnorm(x, pred_i, vcov, logdet), d),
                                         nrow = d,
                                         byrow = TRUE)
                },
                lowerLimit = pred_i - w,
                upperLimit = pred_i + w,
                fDim       = d,
                tol        = rel.acc,
                absError   = 0.0001,
                vectorInterface = TRUE
            )$integral
        })
        
        #Applyign the mask if provided
        if (!is.null(mask)) {vardist[t(mask)] <- NA}
        
        vardist <- apply(vardist, 1, mean, na.rm = TRUE)
        
        vcv <- vcv + diag(vardist)
    }
    
    #Printing the result
    return(vcv)
}

#Computing the Psi vector
QGmvpsi <- function(mu = NULL, 
                    vcov,
                    d.link.inv,
                    predict = NULL,
                    rel.acc = 0.001,
                    width = 10,
                    mask = NULL)
{
    #Setting the integral width according to vcov (lower mean-w, upper mean+w)
    w <- sqrt(diag(vcov)) * width
    
    #Number of dimensions
    d <- length(w)
    
    # Computing the logdet of vcov
    logdet <- calc_logdet(vcov)
    
    #If predict is not included, then use mu, and 
    if(is.null(predict)) { 
        if(is.null(mu)) {
            stop("Please provide either mu or predict.")
        } else {
            predict <- matrix(mu, nrow = 1)
        }
    }
    
    #Computing the mean
    #The double apply is needed to compute the mean for all "predict" values,
    #then average over them
    Psi <- apply(predict, 1, function(pred_i) {
        cubature::hcubature(
            f  = function(x) {
                d.link.inv(x) * matrix(rep(vec_mvnorm(x, pred_i, vcov, logdet), d),
                                       nrow = d,
                                       byrow = TRUE)
            },
            lowerLimit = pred_i - w,
            upperLimit = pred_i + w,
            fDim       = d,
            tol        = rel.acc,
            absError   = 0.0001,
            vectorInterface = TRUE
        )$integral
    })
    
    #Applyign the mask if provided
    if (!is.null(mask)) {Psi[t(mask)] <- NA}
    
    Psi <- apply(Psi, 1, mean, na.rm = TRUE)
    #Make Psi a matrix
    Psi <- diag(Psi)
    #print Psi
    return(Psi)
}

## ----------------- Wrapper function for general calculation ------------------- ##

QGmvparams <- function(mu = NULL,
                       vcv.G,
                       vcv.P,
                       models,
                       predict = NULL,
                       rel.acc = 0.001,
                       width = 10,
                       n.obs = NULL,
                       theta = NULL,
                       verbose = TRUE,
                       mask = NULL)
{
    #Error if ordinal is used (multivariate code not available yet)
    if ("ordinal" %in% models) {
        stop("Multivariate functions of QGglmm are 
             not able to address ordinal traits (yet?).")
    }
    
    #Setting the integral width according to vcov (lower mean-w, upper mean+w)
    w <- sqrt(diag(vcv.P)) * width
    
    #Number of dimensions
    d <- length(w)
    
    #If predict is not included, then use mu, and 
    if(is.null(predict)) {
        if(is.null(mu)) {
            stop("Please provide either mu or predict.")
        } else {
            predict <- matrix(mu, nrow = 1)
        }
    }
    
    #Dimensions checks
    if (length(models) != d | 
        nrow(vcv.G) != d | 
        ncol(vcv.G) != d | 
        nrow(vcv.P) != d | 
        ncol(vcv.P) != d | 
        ncol(predict) != d) 
    {
        stop("Dimensions are incompatible, 
             please check the dimensions of the input.")
    }
    
    #Defining the link/distribution functions
    #If a vector of names were given
    if (!(is.list(models))) {
        if (!is.character(models)) {
            stop("models should be either a list or a vector of characters")
        }
        
        # Need to account for n.obs
        binN <- grepl("binomN", models)
        if (any(binN) & length(n.obs) != sum(binN)) {
            stop("Please provide a number of n.obs equal to the number binomN models.")
        }
        list.n.obs <- rep(list(NULL), length(models))
        list.n.obs[binN] <- n.obs
        
        # Need to account for theta
        negbin <- grepl("negbin", models)
        if (any(negbin) & length(theta) != sum(negbin)) {
            stop("Please provide a number of theta equal to the number binomN models.")
        }
        list.theta <- rep(list(NULL), length(models))
        list.theta[negbin] <- theta
        
        models <- mapply(function(name, n.obs, theta) {
                            QGlink.funcs(name = name, n.obs = n.obs, theta = theta)
                         },
                         name       = models,
                         n.obs      = list.n.obs,
                         theta      = list.theta,
                         SIMPLIFY   = FALSE)
    }
    
    #Now we can compute the needed functions
    inv.links <- function(mat) {
        res <- mat
        for (i in 1:d) {
            res[i, ] <- models[[i]]$inv.link(mat[i, ])
        }
        res
    }
    var.funcs <- function(mat) {
        res <- mat
        for (i in 1:d) {
            res[i, ] <- models[[i]]$var.func(mat[i, ])
        }
        res
    }
    d.inv.links <- function(mat) {
        res <- mat
        for (i in 1:d) {
            res[i, ] <- models[[i]]$d.inv.link(mat[i, ])
        }
        res
    }
    
    #Computing the observed mean
    if (verbose) {
        print("Computing observed mean...")
    }
    z_bar <- QGmvmean(mu        = mu,
                      vcov      = vcv.P,
                      link.inv  = inv.links,
                      predict   = predict,
                      rel.acc   = rel.acc,
                      width     = width,
                      mask      = mask)
    
    #Computing the variance-covariance matrix
    if (verbose) {
        print("Computing variance-covariance matrix...")
    }
    vcv.P.obs <- QGvcov(mu          = mu,
                        vcov        = vcv.P,
                        link.inv    = inv.links,
                        var.func    = var.funcs,
                        mvmean.obs  = z_bar,
                        predict     = predict, 
                        rel.acc     = rel.acc,
                        width       = width,
                        exp.scale   = FALSE,
                        mask        = mask)
    
    if (verbose) {
        print("Computing Psi...")
    }
    
    Psi <- 
        QGmvpsi(mu          = mu,
                vcov        = vcv.P,
                d.link.inv  = d.inv.links,
                predict     = predict, 
                rel.acc     = rel.acc, 
                width       = width,
                mask        = mask)
    
    vcv.G.obs <- Psi %*% vcv.G %*% t(Psi)
    
    #Return a list of QG parameters on the observed scale
    return(list(mean.obs    = z_bar,
                vcv.P.obs   = vcv.P.obs,
                vcv.G.obs   = vcv.G.obs))
}

##  --------------Function to calculate the evolutive prediction-------- ##

QGmvpred <- function(mu = NULL,
                     vcv.G,
                     vcv.P,
                     fit.func,
                     d.fit.func,
                     predict = NULL, 
                     rel.acc = 0.001,
                     width = 10,
                     verbose = TRUE,
                     mask = NULL)
{
    #Setting the integral width according to vcv.P (lower mean-w, upper mean+w)
    w <- sqrt(diag(vcv.P)) * width
    
    #Number of dimensions
    d <- length(w)
    
    # Computing the logdet of vcv.P
    logdet <- calc_logdet(vcv.P)
    
    #If predict is not included, then use mu, and 
    if(is.null(predict)) {
        if(is.null(mu)) {
            stop("Please provide either mu or predict.")
        } else {
            predict <- matrix(mu, nrow = 1)
        }
    }
    
    #Dimensions checks
    if (nrow(vcv.G) != d |
        ncol(vcv.G) != d | 
        nrow(vcv.P) != d | 
        ncol(vcv.P) != d) 
    {
        stop("Dimensions are incompatible, 
             please check the dimensions of the input.")
    }
    
    #Calculating the latent mean fitness
    if (verbose) {
        print("Computing mean fitness...")     
    }
    
    Wbar <- 
        mean(apply(predict, 1,function(pred_i) {
            cubature::hcubature(
                f  = function(x) {
                    fit.func(x) * vec_mvnorm(x, pred_i, vcv.P, logdet)
                },
                lowerLimit = pred_i - w,
                upperLimit = pred_i + w,
                #Note that fDim = 1 because fitness.func yields a scalar
                fDim       = 1,
                tol        = rel.acc,
                absError   = 0.0001,
                vectorInterface = TRUE
            )$integral
        }))
    
    #Calculating the covariance between latent trait and latent fitness
    if (verbose) {
        print("Computing the latent selection and response...")
    }
    
    #Computing the derivative of fitness
    dW <- 
        apply(predict, 1, function(pred_i) {
            cubature::hcubature(
                f  = function(x) {
                    d.fit.func(x) * matrix(rep(vec_mvnorm(x, pred_i, vcv.P, logdet), d),
                                           nrow = d,
                                           byrow = TRUE)
                },
                lowerLimit = pred_i - w,
                upperLimit = pred_i + w,
                fDim       = d,
                tol        = rel.acc,
                absError   = 0.0001,
                vectorInterface = TRUE
            )$integral
        })
    
    #Applyign the mask if provided
    if (!is.null(mask)) {dW[t(mask)] <- NA}
    
    dW <- apply(dW, 1, mean)
    
    #Computing the selection
    if (dim(predict)[1] > 1) {
        sel <- as.vector(((vcv.P + var(predict)) %*% dW) / Wbar) 
    } else {
        sel <- as.vector((vcv.P %*% dW) / Wbar)
    }
    
    #Computing the evolutionary response
    resp <- as.vector((vcv.G %*% dW) / Wbar)
    
    #Returning the results on the latent scale
    return(list(mean.lat.fit = Wbar, 
                lat.grad = dW / Wbar,
                lat.sel = sel,
                lat.resp = resp))
}
