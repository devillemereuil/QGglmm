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

##-------------------------------General functions--------------------------------------------

QGicc <- function(mu=NULL, var.comp, var.p, model="", width=10, predict=NULL, closed.form=TRUE, custom.model=NULL, n.obs=NULL, cut.points=NULL, theta=NULL, verbose=TRUE){
  if(length(mu)>1 | length(var) >1) stop("The parameters mu and var must be of length 1, please check your input.")
  if (is.null(predict)) { if(is.null(mu)) {stop("Please provide either mu or predict.")} else {predict=mu;}}
  
  #Using analytical solutions if possible (and asked for, see closed.form arg)
  if (model=="Gaussian"&closed.form) {
    if (verbose) print("Using the closed forms for a Gaussian model with identity link (e.g. LMM).")
    qg.Gaussian.icc(mu=mu,var.comp=var.comp,var.p=var.p,predict=predict)
  } else if (model=="Poisson.log"&closed.form) {
    if(verbose) print("Using the closed forms for a Poisson-log model.")
    qg.Poisson.log.icc(mu=mu,var.comp=var.comp,var.p=var.p,predict=predict)
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


QGmvicc<-function(mu=NULL,vcov.comp,vcov.p,models,predict=NULL,rel.acc=0.01,width=10,n.obs=NULL,theta=NULL,verbose=TRUE) {
  #Setting the integral width according to vcov (lower mean-w, upper mean+w)
  w1<-sqrt(diag(vcov.p-vcov.comp))*width
  w2<-sqrt(diag(vcov.comp))*width
  #Number of dimensions
  d<-length(w1)
  #If no fixed effects were included in the model
  if (is.null(predict)) predict=matrix(mu,nrow=1)
  
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
  z_bar<-QGmvmean(mu=mu,vcov=vcov.p,link.inv=inv.links,predict=predict,rel.acc=rel.acc,width=width)
  #Computing the variance-covariance matrix
  if (verbose) print("Computing phenotypic variance-covariance matrix...")
  vcv.P.obs<-QGvcov(mu=mu,vcov=vcov.p,link.inv=inv.links,var.func=var.funcs,mvmean.obs=z_bar,predict=predict,rel.acc=rel.acc,width=width,exp.scale=FALSE)
  
  #Function giving the conditional expectancy
  cond_exp <- function(t) {
  apply(#
      apply(predict,1,
            function(pred_i){
	      cuhre(ndim=d,ncomp=d,
              integrand=function(x){inv.links(x)*dmvnorm(x,pred_i+t,vcov.p-vcov.comp)},
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
	  (cond_exp(x)%*%t(cond_exp(x)))[upper.tri(x%*%t(x),diag=TRUE)]*dmvnorm(x,rep(0,d),vcov.comp)
     },
     lower=-w2,upper=w2,rel.tol=rel.acc,abs.tol= 0.001,
     flags=list(verbose=0))$value
  
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
