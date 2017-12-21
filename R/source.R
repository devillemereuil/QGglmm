#----------------------------------------------------------------
#	QGglmm R package
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

##-------------------------- General functions ---------------------- ##



#Calculating the observed/expected scale mean
QGmean <- function(mu = NULL, 
                   var, 
                   link.inv, 
                   predict = NULL, 
                   width = 10) 
{
    if(length(mu) > 1 | length(var) > 1) {
        stop("The parameters mu and var must be of length 1, 
             please check your input.")
    }
    #If no fixed effects were included in the model
    if (is.null(predict)) {
        if(is.null(mu)) {
            stop("Please provide either mu or predict.")
        } else {
            predict <- mu
        }
    }
    mean(
        sapply(predict,
               function(pred_i) {
                   integrate(f = function(x) {
                                link.inv(x) * dnorm(x, pred_i, sqrt(var))
                             },
                             lower = pred_i - width * sqrt(var),
                             upper = pred_i + width * sqrt(var)
                            )$value
                   }
               )
        )
}

#Calculating the expected scale variance
QGvar.exp <- function(mu = NULL,
                      var,
                      link.inv,
                      obs.mean = NULL,
                      predict = NULL,
                      width = 10) 
{
    if(length(mu) > 1 | length(var) > 1) {
        stop("The parameters mu and var must be of length 1, 
             please check your input.")
    }
    #If no fixed effects were included in the model
    if (is.null(predict)) {
        if(is.null(mu)) {
            stop("Please provide either mu or predict.")
        } else {
            predict <- mu
        }
    }
    #If not provided, compute the obsereved mean
    if (is.null(obs.mean)) {
        obs.mean <- QGmean(mu        = mu, 
                           var       = var, 
                           link.inv  = link.inv,
                           width     = width,
                           predict   = predict)
    }
    #Note: Using Koenig's formula
    mean(
        sapply(predict,
               function(pred_i) {
                   integrate(f = function(x) {
                                    (link.inv(x)^2) * 
                                        dnorm(x, pred_i, sqrt(var))
                                 }, 
                             lower = pred_i - width * sqrt(var), 
                             upper = pred_i + width * sqrt(var)
                   )$value
               }
        )
    ) - (obs.mean^2)
}

#Calculating the "distribution" variance
QGvar.dist <- function(mu = NULL, 
                       var, 
                       var.func, 
                       predict = NULL, 
                       width = 10) 
{
    if(length(mu) > 1 | length(var) > 1) {
        stop("The parameters mu and var must be of length 1, 
             please check your input.")
    }
    #If no fixed effects were included in the model
    if (is.null(predict)) {
        if(is.null(mu)) {
            stop("Please provide either mu or predict.")
        } else {
            predict <- mu
        }
    }
    mean(
        sapply(predict, 
               function(pred_i) {
                   integrate(f = function(x) {
                                    var.func(x) * dnorm(x, pred_i, sqrt(var))
                                 },
                             lower = pred_i - width * sqrt(var), 
                             upper = pred_i + width * sqrt(var)
                   )$value
               }
        )
    )
}

#Calculating "psi" for the observed additive genetic variance computation
QGpsi <- function(mu = NULL, 
                  var, 
                  d.link.inv, 
                  predict = NULL, 
                  width = 10)
{
    if(length(mu) > 1 | length(var) > 1) {
        stop("The parameters mu and var must be of length 1, 
             please check your input.")
    }
    #If no fixed effects were included in the model
    if (is.null(predict)) {
        if(is.null(mu)) {
            stop("Please provide either mu or predict.")
        } else {
            predict <- mu
        }
    }
    mean(
        sapply(predict, 
               function(pred_i) {
                   integrate(f = function(x) {
                                    d.link.inv(x) * dnorm(x, pred_i, sqrt(var))
                                 },
                             lower = pred_i - width * sqrt(var), 
                             upper = pred_i + width * sqrt(var)
                   )$value
               }
        )
    )
}

##  ------------------ Dictionary functions ----------------- ##

#Function creating the needed functions according to the link name
QGlink.funcs <- function(name, 
                         n.obs = NULL,
                         theta = NULL)
{
    if (name == "Gaussian") {
        inv.link    <- function(x) {x}
        var.func    <- function(x) {0}
        d.inv.link  <- function(x) {1}
    } else if (name == "binom1.probit") {
        inv.link    <- function(x) {pnorm(x)}
        var.func    <- function(x) {pnorm(x) * (1 - pnorm(x))}
        d.inv.link  <- function(x) {dnorm(x)}
    } else  if (name == "binomN.probit") {
        if (is.null(n.obs)) {
            stop("binomN model used, 
                 but no observation number (n.obs) defined.")
        }
        inv.link    <- function(x) {n.obs * pnorm(x)}
        var.func    <- function(x) {n.obs * pnorm(x) * (1 - pnorm(x))}
        d.inv.link  <- function(x) {n.obs * dnorm(x)}
    } else if (name == "binom1.logit") {
        inv.link    <- function(x) {plogis(x)}
        var.func    <- function(x) {plogis(x) * (1 - plogis(x))}
        d.inv.link  <- function(x) {dlogis(x)}
    } else if (name == "binomN.logit") {
        if (is.null(n.obs)) {
            stop("binomN model used, 
                 but no observation number (n.obs) defined.")
        }
        inv.link    <- function(x) {n.obs * plogis(x)}
        var.func    <- function(x) {n.obs * plogis(x) * (1 - plogis(x))}
        d.inv.link  <- function(x) {n.obs * dlogis(x)}
    } else if (name == "ordinal") {
        stop("ordinal models are particular, 
             please use the QGparams function only")
    } else if (name == "Poisson.log") {
        inv.link    <- function(x) {exp(x)}
        var.func    <- function(x) {exp(x)}
        d.inv.link  <- function(x) {exp(x)}
    } else if (name == "Poisson.sqrt") {
        inv.link    <- function(x) {x^2}
        var.func    <- function(x) {x^2}
        d.inv.link  <- function(x) {2 * x}
    } else if (name == "negbin.log") {
        if (is.null(theta)) {stop("negbin model used, but theta not defined.")}
        inv.link    <- function(x) {exp(x)}
        var.func    <- function(x) {exp(x) + (exp(2 * x) / theta)}
        d.inv.link  <- function(x) {exp(x)}
    } else if (name == "negbin.sqrt") {
        if (is.null(theta)) {
            stop("negbin model used, but theta not defined.")
        }
        inv.link    <- function(x) {x^2}
        var.func    <- function(x) {(x^2) + ((x^4) / theta)}
        d.inv.link  <- function(x) {2 * x}
    } else {
        stop("Invalid model name. 
             Use a valid model name or enter a custom model specification.")
    }
    list(inv.link = inv.link, var.func = var.func, d.inv.link = d.inv.link)
}

## -----------Special functions for known analytical solutions-------------- ##

qg.Gaussian <- function(mu = NULL, 
                        var.a, 
                        var.p, 
                        predict = NULL) 
{
    if(length(mu) > 1 | length(var.a) != 1 | length(var.p) != 1) {
        stop("The parameters mu and var must be of length 1, 
             please check your input.")
    }
    #If no fixed effects were included in the model
    if (is.null(predict)) {
        if(is.null(mu)) {
            stop("Please provide either mu or predict.")
        } else {
            predict <- mu
        }
    }
    
    #Nothing to be done, except averaging over predict
    if (length(predict) == 1) {
        var_fixed <- 0
    } else {
        var_fixed <- var(predict)
    }
    
    data.frame(mean.obs     = mean(predict),
               var.obs      = var.p + var_fixed, 
               var.a.obs    = var.a, 
               h2.obs       = var.a / (var.p + var_fixed))
}

qg.binom1.probit <- function(mu = NULL, 
                             var.a, 
                             var.p, 
                             predict = NULL) 
{
    if(length(mu) > 1 | length(var.a) != 1 | length(var.p) != 1) {
        stop("The parameters mu and var must be of length 1, 
             please check your input.")
    }
    #If no fixed effects were included in the model
    if (is.null(predict)) {
        if(is.null(mu)) {
            stop("Please provide either mu or predict.")
        } else {
            predict <- mu
        }
    }
    
    #Observed mean
    p <- mean(1 - pnorm(0, predict, sqrt(var.p + 1)))
    
    #Observed variance
    var_obs <- p * (1 - p)
    
    #Psi
    Psi <- mean(dnorm(0, (predict), sqrt(var.p + 1)))
    
    data.frame(mean.obs     = p,
               var.obs      = var_obs, 
               var.a.obs    = (Psi^2) * var.a, 
               h2.obs       = ((Psi^2) * var.a) / var_obs)
}

qg.binomN.probit <- function(mu = NULL,
                             var.a,
                             var.p,
                             n.obs,
                             predict = NULL,
                             width = 10) 
{
    if(length(mu) > 1 | length(var.a) != 1 | length(var.p) != 1) {
        stop("The parameters mu and var must be of length 1, 
             please check your input.")
    }
    #If no fixed effects were included in the model
    if (is.null(predict)) {
        if(is.null(mu)) {
            stop("Please provide either mu or predict.")
        } else {
            predict <- mu
        }
    }
  
  #Observed mean
  p <- n.obs * mean(1 - pnorm(0, predict, sqrt(var.p + 1)))
  
  #Observed variance
  prob.sq.int <- 
      mean(
          sapply(predict, 
                 function(pred_i) {
                     integrate(f = function(x) {
                                       (pnorm(x)^2) * 
                                           dnorm(x, pred_i, sqrt(var.p))
                                    },
                               lower = pred_i - width * sqrt(var.p),
                               upper = pred_i + width * sqrt(var.p)
                     )$value
                 }
          )
      )
  var_obs <- ((n.obs^2) - n.obs) * prob.sq.int - p^2 +p
  
  #Psi
  Psi <- n.obs * mean(dnorm(0, (predict), sqrt(var.p + 1)))
  
  data.frame(mean.obs   = p, 
             var.obs    = var_obs, 
             var.a.obs  = (Psi^2) * var.a, 
             h2.obs     = ((Psi^2) * var.a) / var_obs)
}

qg.Poisson.log <- function(mu = NULL, 
                           var.a, 
                           var.p, 
                           predict = NULL) 
{
    if(length(mu) > 1 | length(var.a) != 1 | length(var.p) != 1) {
        stop("The parameters mu and var must be of length 1, 
             please check your input.")
    }
    #If no fixed effects were included in the model
    if (is.null(predict)) {
        if(is.null(mu)) {
            stop("Please provide either mu or predict.")
        } else {
            predict <- mu
        }
    }
    
    #Observed mean
    lambda <- mean(exp(predict + (var.p / 2)))
    
    #Mean of lambda square, needed for the following
    lambda_sq <- mean(exp(2 * (predict + var.p / 2)))
    
    #Observed variance
    var_obs <- lambda_sq * exp(var.p) - lambda^2 + lambda
    
    data.frame(mean.obs     = lambda, 
               var.obs      = var_obs, 
               var.a.obs    = (lambda^2) * var.a,
               h2.obs       = ((lambda^2) * var.a) / var_obs)
}

qg.Poisson.sqrt <- function(mu = NULL,
                            var.a,
                            var.p,
                            predict = NULL)
{
    if(length(mu) > 1 | length(var.a) != 1 | length(var.p) != 1) {
        stop("The parameters mu and var must be of length 1, 
             please check your input.")
    }
    #If no fixed effects were included in the model
    if (is.null(predict)) {
        if(is.null(mu)) {
            stop("Please provide either mu or predict.")
        } else {
            predict <- mu
        }
    }
    
    #Observed mean
    lambda <- mean((predict)^2 + var.p)
    
    #Observed variance
    var_obs <- mean((predict)^4 + 
                    6 * var.p * ((predict)^2) + 
                    3 * (var.p^2)) - 
               lambda^2 + lambda
    
    #Psi
    Psi <- mean(2 * (predict))
    
    data.frame(mean.obs     = lambda, 
               var.obs      = var_obs, 
               var.a.obs    = (Psi^2) * var.a, 
               h2.obs       = ((Psi^2) * var.a) / var_obs)
}

qg.negbin.log <- function(mu = NULL,
                          var.a,
                          var.p, 
                          theta, 
                          predict = NULL) 
{
    if(length(mu) > 1 | length(var.a) != 1 | length(var.p) != 1) {
        stop("The parameters mu and var must be of length 1, 
             please check your input.")
    }
    #If no fixed effects were included in the model
    if (is.null(predict)) {
        if(is.null(mu)) {
            stop("Please provide either mu or predict.")
        } else {
            predict <- mu
        }
    }
    
    #Observed mean
    lambda <- mean(exp(predict + (var.p / 2)))
    
    #Mean of lambda square, needed for the following
    lambda_sq <- mean(exp(2 * (predict + var.p / 2)))
    
    #Observed variance
    var_obs <- lambda_sq * exp(var.p) - 
               lambda^2 + 
               lambda + 
               mean(exp(2 * (predict + var.p))) / theta
    
    data.frame(mean.obs     = lambda,
               var.obs      = var_obs,
               var.a.obs    = (lambda^2) * var.a, 
               h2.obs       = ((lambda^2) * var.a) / var_obs)
}

qg.negbin.sqrt <- function(mu = NULL,
                           var.a,
                           var.p,
                           theta,
                           predict = NULL)
{
    if(length(mu) > 1 | length(var.a) != 1 | length(var.p) != 1) {
        stop("The parameters mu and var must be of length 1, 
             please check your input.")
    }
    #If no fixed effects were included in the model
    if (is.null(predict)) {
        if(is.null(mu)) {
            stop("Please provide either mu or predict.")
        } else {
            predict <- mu
        }
    }
    
    #Observed mean
    lambda <- mean((predict)^2 + var.p)
    
    #Observed variance
    var_obs <- mean((predict)^4 + 
                    6 * var.p * ((predict)^2) + 
                    3 * (var.p^2)) -
               lambda^2 + lambda + 
               (mean((predict)^4 + 
                     6 * var.p * ((predict)^2) + 
                     3 * (var.p^2)) / theta)
    
    #Psi
    Psi <- mean(2 * (predict))
    
    data.frame(mean.obs     = lambda, 
               var.obs      = var_obs, 
               var.a.obs    = (Psi^2) * var.a, 
               h2.obs       = ((Psi^2) * var.a) / var_obs)
}

##------------------Meta-function for general calculation------------------- ##

QGparams <- function(mu = NULL, 
                     var.a, 
                     var.p, 
                     model = "", 
                     width = 10, 
                     predict = NULL, 
                     closed.form = TRUE, 
                     custom.model = NULL, 
                     n.obs = NULL, 
                     cut.points = NULL, 
                     theta = NULL, 
                     verbose = TRUE)
{
    if(length(mu) > 1 | length(var.a) != 1 | length(var.p) != 1) {
        stop("The parameters mu and var must be of length 1, 
             please check your input.")
    }
    #If no fixed effects were included in the model
    if (is.null(predict)) {
        if(is.null(mu)) {
            stop("Please provide either mu or predict.")
        } else {
            predict <- mu
        }
    }
    
    ## Using analytical solutions if possible (and asked for: see closed.form)
    if (model == "Gaussian" & closed.form) {   
        #Gaussian model with identity link (e.g. LMM)
        if (verbose) {
            print("Using the closed forms for a Gaussian model
                  with identity link (e.g. LMM).")
        }
        
        qg.Gaussian(mu      = mu, 
                    var.a   = var.a, 
                    var.p   = var.p, 
                    predict = predict)
        
    } else if (model == "binom1.probit" & closed.form) {
        #Binary.probit model
        if (verbose) {
            print("Using the closed forms for a Binomial1 - probit model.")
        }
        
        qg.binom1.probit(mu         = mu,
                         var.a      = var.a,
                         var.p      = var.p,
                         predict    = predict)
        
    } else if (model == "threshold") {
        #Binary.probit model
        if (verbose) {
            print("Using the closed forms for a threshold model 
                  (e.g. for MCMCglmm package). 
                  Closed form argument is ignored...")
        }
        
        qg.binom1.probit(mu         = mu,
                         var.a      = var.a,
                         var.p      = var.p - 1,
                         predict    = predict)
        
    } else if (model == "binomN.probit" & closed.form) {
        #Binomial - not - binary model
        if (is.null(n.obs)) {
            stop("binomN.probit model used, 
                 but no observation number (n.obs) defined.")
        }
        if (verbose) {
            print("Using a semi - closed form for a BinomialN - probit model.")
        }
        
        warning("Some component (var.obs) are not totally 
                from a closed form solution: an integral is computed")
        
        qg.binomN.probit(mu         = mu, 
                         var.a      = var.a, 
                         var.p      = var.p, 
                         predict    = predict,
                         n.obs      = n.obs, 
                         width      = width)
        
    } else if (model == "Poisson.log" & closed.form){
        #Poisson - log model
        if(verbose) {
            print("Using the closed forms for a Poisson - log model.")
        }
        
        qg.Poisson.log(mu       = mu,
                       var.a    = var.a,
                       var.p    = var.p,
                       predict  = predict)
        
    } else if (model == "Poisson.sqrt" & closed.form){
        #Poisson - sqrt model
        if(verbose) {
            print("Using the closed forms for a Poisson - sqrt model.")
        }
        
        qg.Poisson.sqrt(mu      = mu,
                        var.a   = var.a,
                        var.p   = var.p,
                        predict = predict)
        
    } else if (model == "negbin.log" & closed.form){
        #NegBin - log model
        if (is.null(theta)) {
            stop("negbin model used, but theta not defined.")
        }
        if(verbose) {
            print("Using the closed forms for a NegativeBinomial - log model.")
        }
        
        qg.negbin.log(mu        = mu, 
                      var.a     = var.a, 
                      var.p     = var.p, 
                      predict   = predict, 
                      theta     = theta)
        
    } else if (model == "negbin.sqrt" & closed.form){
        #NegBin - sqrt model
        if (is.null(theta)) {
            stop("negbin model used, but theta not defined.")
        }
        if(verbose) {
            print("Using the closed forms for a NegativeBinomial - sqrt model.")
        }
        
        qg.negbin.sqrt(mu       = mu,
                       var.a    = var.a,
                       var.p    = var.p, 
                       predict  = predict, 
                       theta    = theta)
        
    } else if (model == "ordinal"){
        #Ordinal model
        if (is.null(cut.points)) {
            stop("cut points must be specified to use the ordinal model.")
        }
        if(verbose) {
            print("Using the closed forms for an ordinal model 
                  (ignoring the closed.form argument)")
        }
        
        qg.ordinal(mu       = mu,
                   var.a    = var.a, 
                   var.p    = var.p,
                   predict  = predict, 
                   cut.points = cut.points)
    } else {
        
        ##Else, use the general integral equations
        
        #Use a custom model if defined, otherwise look into the "Dictionary"
        if (is.null(custom.model)) {
            if (model == "") {
                stop("The function requires either model or custom.model.")
            } else {
                funcs <- QGlink.funcs(model, n.obs = n.obs, theta = theta)
            }
        } else {
            funcs <- custom.model
        }
        
        #Observed mean computation
        if (verbose) {
            print("Computing observed mean...")
        }
        z_bar <- QGmean(mu          = mu,
                        var         = var.p, 
                        link.inv    = funcs$inv.link,
                        width       = width,
                        predict     = predict)
        
        #Variances computation
        if (verbose) {
            print("Computing variances...")
        }
        var_exp <- QGvar.exp(mu         = mu,
                             var        = var.p,
                             link.inv   = funcs$inv.link,
                             obs.mean   = z_bar, 
                             width      = width,
                             predict    = predict)
        var_dist <- QGvar.dist(mu       = mu, 
                               var      = var.p, 
                               var.func = funcs$var.func,
                               width    = width,
                               predict  = predict)
        var_obs <- var_exp + var_dist
        
        #Psi computation (for the observed additive genetic variance)
        if (verbose) {
            print("Computing Psi...")
        }
        Psi <- QGpsi(mu         = mu,
                     var        = var.p, 
                     d.link.inv = funcs$d.inv.link,
                     width      = width,
                     predict    = predict)
        
        #Return a data.frame with the calculated components
        data.frame(mean.obs     = z_bar,
                   var.obs      = var_obs,
                   var.a.obs    = (Psi^2) * var.a,
                   h2.obs       = ((Psi^2) * var.a) / var_obs)
}
}

## -----------Function to calculate the evolutive prediction--------------- ##

QGpred <- function(mu = NULL, 
                   var.a, 
                   var.p, 
                   fit.func, 
                   d.fit.func, 
                   width = 10, 
                   predict = NULL, 
                   verbose = TRUE) 
{
    if(length(mu) > 1 | length(var.a) != 1 | length(var.p) != 1) {
        stop("The parameters mu, var.a and var.p must be of length 1, 
             please check your input.")
    }
    if (is.null(predict)) { 
        if(is.null(mu)) {
            stop("Please provide either mu or predict.")
        } else {
            predict <- mu
        }
    }
    
    #Calculating the latent mean fitness
    if (verbose) {
        print("Computing mean fitness...")
    }
    
    Wbar <- 
        mean(
            sapply(predict, 
                   function(pred_i) {
                       integrate(f = function(x) {
                                        fit.func(x) * 
                                            dnorm(x, pred_i, sqrt(var.p))
                                     },
                                 lower = pred_i - width * sqrt(var.p),
                                 upper = pred_i + width * sqrt(var.p)
                       )$value
                   }
            )
        )
    
    #Calculating the covariance between latent trait and latent fitness
    if (verbose) {
        print("Computing the latent selection and response...")
    }
    
    #Computing the derivative of fitness
    dW <- 
        mean(
            sapply(predict, 
                   function(pred_i) {
                       integrate(f = function(x) {
                                        d.fit.func(x) * 
                                            dnorm(x, pred_i, sqrt(var.p))
                                     },
                                 lower = pred_i - width * sqrt(var.p),
                                 upper = pred_i + width * sqrt(var.p)
                       )$value
                   }
            )
        )
    
    #Computing the selection
    if (length(predict) > 1) {
        sel <- (var.p + var(predict)) * dW / Wbar
    } else {
        sel <- var.p * dW / Wbar
    } 
    
    #Computing the evolutionary response
    resp <- var.a * dW / Wbar
    
    #Returning the results on the latent scale
    data.frame(mean.lat.fit = Wbar, 
               lat.grad     = dW / Wbar,
               lat.sel      = sel,
               lat.resp     = resp)
}
