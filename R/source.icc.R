#----------------------------------------------------------------
#	QGglmm R package (functions for ICC estimates)
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
#    GNU General Public License for more detailsls.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

##-----------------------------Special functions for known analytical solutions--------------

qg.Gaussian.icc=function(mu=NULL,var.comp,var.p,predict=NULL) {
  if(length(mu)>1 | length(var.comp)!=1 | length(var.p) != 1) stop("The parameters mu, var.comp and var.p must be of length 1, please check your input.")
  if (is.null(predict)) { if(is.null(mu)) {stop("Please provide either mu or predict.")} else {predict=mu;}}
  #Nothing to be done, except averaging over predict
  if (length(predict)==1) {var_fixed=0} else {var_fixed=var(predict)}
  data.frame(mean.obs=mean(predict),var.obs=var.p+var_fixed,var.comp.obs=var.comp,icc.obs=var.comp/(var.p+var_fixed))
}

qg.binom1.probit.icc <- function(mu=NULL,var.comp,var.p,predict=NULL, width=10){
  if(length(mu)>1 | length(var.comp)!=1 | length(var.p) != 1) stop("The parameters mu, var.comp and var.p must be of length 1, please check your input.")
  if (is.null(predict)) { if(is.null(mu)) {stop("Please provide either mu or predict.")} else {predict=mu;}}
  #Observed mean
  p=mean(1-pnorm(0,predict,sqrt(var.p+1)))
  #Observed variance
  var_obs=p*(1-p)
  #Component variance
  var_comp_obs = integrate(f=function(x){
    sapply(x,function(x_i){
      (mean(pnorm(x_i,-predict,sqrt(var.p-var.comp+1)))^2)*dnorm(x_i,0,sqrt(var.comp))
    })
  },lower=-width*sqrt(var.comp),upper=width*sqrt(var.comp))$value - (p^2)
  data.frame(mean.obs=p,var.obs=var_obs,var.comp.obs=var_comp_obs,icc.obs=var_comp_obs/var_obs)
}

qg.binomN.probit.icc <- function(mu=NULL,var.comp,var.p,n.obs,predict=NULL, width=10){
  if(length(mu)>1 | length(var.comp)!=1 | length(var.p) != 1) stop("The parameters mu, var.comp and var.p must be of length 1, please check your input.")
  if (is.null(predict)) { if(is.null(mu)) {stop("Please provide either mu or predict.")} else {predict=mu;}}
  #Observed mean
  p=n.obs*mean(1-pnorm(0,predict,sqrt(var.p+1)))
  #Observed variance
  prob.sq.int=mean(sapply(predict,function(pred_i){integrate(f=function(x){(pnorm(x)**2)*dnorm(x,pred_i,sqrt(var.p))},lower=pred_i-width*sqrt(var.p),upper=pred_i+width*sqrt(var.p))$value}))
  var_obs=((n.obs**2)-n.obs)*prob.sq.int - p**2 +p
  #Component variance
  var_comp_obs = integrate(f=function(x){
    sapply(x,function(x_i){
      (mean(n.obs*pnorm(x_i,-predict,sqrt(var.p-var.comp+1)))^2)*dnorm(x_i,0,sqrt(var.comp))
    })
  },lower=-width*sqrt(var.comp),upper=width*sqrt(var.comp))$value - (p^2)
  data.frame(mean.obs=p,var.obs=var_obs,var.comp.obs=var_comp_obs,icc.obs=var_comp_obs/var_obs)
}

qg.Poisson.log.icc <- function(mu=NULL,var.comp,var.p,predict=NULL){
  if(length(mu)>1 | length(var.comp)!=1 | length(var.p) != 1) stop("The parameters mu, var.comp and var.p must be of length 1, please check your input.")
  if (is.null(predict)) { if(is.null(mu)) {stop("Please provide either mu or predict.")} else {predict=mu;}}
  #Observed mean
  lambda = mean(exp(predict+(var.p/2)))
  #Mean of lambda square, needed for the following
  lambda_sq = mean(exp(2*(predict+var.p/2)))
  #Observed variance
  var_obs = lambda_sq*exp(var.p)-lambda**2+lambda
  #Component variance
  var_comp_obs = exp(var.comp + var.p)*(mean(exp(predict))^2) - (lambda^2)
  data.frame(mean.obs=lambda,var.obs=var_obs,var.comp.obs=var_comp_obs,icc.obs=var_comp_obs/var_obs)
}

qg.negbin.log.icc <- function(mu=NULL,var.comp,var.p,theta,predict=NULL){
  if(length(mu)>1 | length(var.comp)!=1 | length(var.p) != 1) stop("The parameters mu, var.comp and var.p must be of length 1, please check your input.")
  if (is.null(predict)) { if(is.null(mu)) {stop("Please provide either mu or predict.")} else {predict=mu;}}
  #Observed mean
  lambda=mean(exp(predict+(var.p/2)))
  #Mean of lambda square, needed for the following
  lambda_sq=mean(exp(2*(predict+var.p/2)))
  #Observed variance
  var_obs=lambda_sq*exp(var.p)-lambda**2+lambda+(mean(exp(2*(predict+var.p)))/theta)
  #Component variance
  var_comp_obs = exp(var.comp + var.p)*(mean(exp(predict))^2) - (lambda^2)
  data.frame(mean.obs=lambda,var.obs=var_obs,var.comp.obs=var_comp_obs,icc.obs=var_comp_obs/var_obs)
}

##-------------------------------General functions--------------------------------------------

QGicc <- function(mu=NULL, var.comp, var.p, model="", width=10, predict=NULL, closed.form=TRUE, custom.model=NULL, n.obs=NULL, theta=NULL, verbose=TRUE){
  #Error if ordinal is used (multivariate code not available yet)
  if ("ordinal" %in% model) {stop("ICC computations are not able to address ordinal traits (yet?).")}
  if(length(mu)>1 | length(var) >1) stop("The parameters mu and var must be of length 1, please check your input.")
  if (is.null(predict)) { if(is.null(mu)) {stop("Please provide either mu or predict.")} else {predict=mu;}}
  
  #Using analytical solutions if possible (and asked for, see closed.form arg)
  if (model=="Gaussian"&closed.form) {
    if (verbose) print("Using the closed forms for a Gaussian model with identity link (e.g. LMM).")
    qg.Gaussian.icc(mu=mu,var.comp=var.comp,var.p=var.p,predict=predict)
  } else if (model=="binom1.probit"&closed.form) {
    if(verbose) print("Using a semi-closed form for a binom1.probit model.")
    warning("Some component (var.comp.obs) are not totally from a closed form solution: an integral is computed")
    qg.binom1.probit.icc(mu=mu,var.comp=var.comp,var.p=var.p,predict=predict,width=width)
  } else if (model=="binomN.probit"&closed.form) {					#Binomial-not-binary model
      if (is.null(n.obs)) {stop("binomN.probit model used, but no observation number (n.obs) defined.")}
      if (verbose) print("Using a semi-closed form for a BinomialN-probit model.")
      warning("Some component (var.obs and var.comp.obs) are not totally from a closed form solution: an integral is computed")
      qg.binomN.probit.icc(mu=mu,var.comp=var.comp,var.p=var.p,predict=predict,n.obs=n.obs,width=width)
  } else if (model=="Poisson.log"&closed.form) {
    if(verbose) print("Using the closed forms for a Poisson-log model.")
    qg.Poisson.log.icc(mu=mu,var.comp=var.comp,var.p=var.p,predict=predict)
  } else if (model=="negbin.log"&closed.form){						#NegBin-log model
      if (is.null(theta)) {stop("negbin model used, but theta not defined.")}
      if(verbose) print("Using the closed forms for a NegativeBinomial-log model.")
      qg.negbin.log.icc(mu=mu,var.comp=var.comp,var.p=var.p,predict=predict,theta=theta)
  } else {
  
    #Use a custom model if defined, otherwise look into the "Dictionary"
    if (is.null(custom.model)) {
      if (model==""){stop("The function requires either model or custom.model.")} else {
	funcs=QGlink.funcs(model,n.obs=n.obs,theta=theta)
      }} else {funcs=custom.model}

    #Observed mean computation
    if (verbose) print("Computing observed mean...")
    z_bar=QGmean(mu=mu,var=var.p,link.inv=funcs$inv.link,width=width,predict=predict)

    #Phenotypic variances computation
    if (verbose) print("Computing phenotypic variance...")
    var_exp=QGvar.exp(mu=mu,var=var.p,link.inv=funcs$inv.link,obs.mean=z_bar,width=width,predict=predict)
    var_dist=QGvar.dist(mu=mu,var=var.p,var.func=funcs$var.func,width=width,predict=predict)
    var_obs=var_exp+var_dist
    
    #Computing conditional variance (i.e. without var.comp)
    var.cond <- var.p - var.comp
    
    #Function giving the conditional expectancy
    cond_exp <- function(t){
      sapply(t,function(t_i){
	mean(sapply(predict,function(pred_i){
	  integrate(f=function(u){funcs$inv.link(u)*dnorm(u,pred_i+t_i,sqrt(var.cond))},
	  lower=pred_i+t_i-width*sqrt(var.cond),upper=pred_i+t_i+width*sqrt(var.cond))$value
	}))
      })
    }
    
    #Computing variance of component on the observed data scale
    #Note: Using Koenig's formula
    if (verbose) print("Computing component variance...")
    var_comp_obs <- integrate(f=function(x){((cond_exp(x)-z_bar)^2)*dnorm(x,0,sqrt(var.comp))},lower=-width*sqrt(var.comp),upper=width*sqrt(var.comp))$value
    
    #Return a data.frame with the calculated components
    data.frame(mean.obs=z_bar,var.obs=var_obs,var.comp.obs=var_comp_obs,icc.obs=var_comp_obs/var_obs)
  }
}


QGmvicc<-function(mu=NULL,vcv.comp,vcv.P,models,predict=NULL,rel.acc=0.01,width=10,n.obs=NULL,theta=NULL,verbose=TRUE,mask=NULL) {
  #Error if ordinal is used (multivariate code not available yet)
  if ("ordinal" %in% models) {stop("Multivariate functions of QGglmm are not able to address ordinal traits (yet?).")}
  #Setting the integral width according to vcov (lower mean-w, upper mean+w)
  w1<-sqrt(diag(vcv.P-vcv.comp))*width
  w2<-sqrt(diag(vcv.comp))*width
  #Number of dimensions
  d<-length(w1)
  #If no fixed effects were included in the model
  if (is.null(predict)) predict<-matrix(mu,nrow=1)
  
  #Defining the link/distribution functions
  #If a vector of names were given
  if (!(is.list(models))) {
    tmp<-models
    models<-lapply(tmp,function(string){QGlink.funcs(string)})
  }
  #Now we can compute the needed functions
  inv.links=function(x){
    res=numeric(d)
    for (i in 1:d) {
      res[i] = models[[i]]$inv.link(x[i])
    }
    res
  }
  var.funcs=function(x){
    res=numeric(d)
    for (i in 1:d) {
      res[i] = models[[i]]$var.func(x[i])
    }
    res
  }
  d.inv.links=function(x){
    res=numeric(d)
    for (i in 1:d) {
      res[i] = models[[i]]$d.inv.link(x[i])
    }
    res
  }
  
  #Computing the observed mean
  if (verbose) print("Computing observed mean...")
  z_bar<-QGmvmean(mu=mu,vcov=vcv.P,link.inv=inv.links,predict=predict,rel.acc=rel.acc,width=width,mask=mask)
  #Computing the variance-covariance matrix
  if (verbose) print("Computing phenotypic variance-covariance matrix...")
  vcv.P.obs<-QGvcov(mu=mu,vcov=vcv.P,link.inv=inv.links,var.func=var.funcs,mvmean.obs=z_bar,predict=predict,rel.acc=rel.acc,width=width,exp.scale=FALSE,mask=mask)
  
  #Function giving the conditional expectancy
  cond_exp <- function(t) {
  apply(#
      apply(predict,1,
            function(pred_i){
	      cuhre(ndim=d,ncomp=d,
              integrand=function(x){inv.links(x)*dmvnorm(x,pred_i+t,vcv.P-vcv.comp)},
              lower=pred_i+t-w1,upper=pred_i+t+w1,rel.tol=rel.acc,abs.tol= 0.001,
              flags=list(verbose=0))$value
	    }
      ),
     1,mean)
  }
  
  if (verbose) print("Computing component variance-covariance matrix...")
  #Computing the upper-triangle matrix of "expectancy of the square"
  v<-cuhre(ndim=d,ncomp=(d^2+d)/2,
      integrand=function(x){
	  (cond_exp(x)%*%t(cond_exp(x)))[upper.tri(x%*%t(x),diag=TRUE)]*dmvnorm(x,rep(0,d),vcv.comp)
     },
     lower=-w2,upper=w2,rel.tol=rel.acc,abs.tol= 0.001,
     flags=list(verbose=0))$value
  #Applying the mask if provided
  if(!is.null(mask)) {
      mask2 <- matrix(FALSE, ncol = ncol(v), nrow = nrow(v))
      mask2[((1:d) * ((1:d) + 1)) / 2, ] <- t(mask)
      v[mask2] <- NA
  }
  
  #Creating the VCV matrix
  vcv<-matrix(NA,d,d)
  vcv[upper.tri(vcv,diag=TRUE)]<-v
  vcv[lower.tri(vcv)]<-vcv[upper.tri(vcv)]

  #Computing the VCV matrix using Keonig's formuka
  vcv <- vcv - z_bar%*%t(z_bar)
  
  #Printing the result
  vcv
  
  #Return a list of QG parameters on the observed scale
  list(mean.obs=z_bar,vcv.P.obs=vcv.P.obs,vcv.comp.obs=vcv)
} 
