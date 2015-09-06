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

##---------------------------------General functions----------------------------------------

#Calculating the observed/expected scale mean (multivariate)
QGmvmean<-function(mu,vcov,link.inv,predict=NULL,rel.acc=0.01,width=10) {
  #Setting the integral width according to vcov (lower mean-w, upper mean+w)
  w<-sqrt(diag(vcov))*width
  #Number of dimensions
  d<-length(w)
  #If predict is not included, then use mu, and 
  if(is.null(predict)) predict=matrix(mu,nrow=1)
  #Computing the mean
  #The double apply is needed to compute the mean for all "predict" values,
  #then average over them
  apply(apply(predict,1,function(pred_i){
      cuhre(ndim=d,ncomp=d,
            integrand=function(x){link.inv(x)*dmvnorm(x,pred_i,vcov)},
            lower=pred_i-w,upper=pred_i+w,rel.tol=rel.acc,abs.tol=0.0001,
            flags=list(verbose=0))$value}),1,mean)
}

#Calculating the expected scale variance-covariance matrix
QGvcov<-function(mu,vcov,link.inv,var.func,mvmean.obs=NULL,predict=NULL,rel.acc=0.01,width=10,exp.scale=FALSE) {
  #Setting the integral width according to vcov (lower mean-w, upper mean+w)
  w<-sqrt(diag(vcov))*width
  #Number of dimensions
  d<-length(w)
  #If no fixed effects were included in the model
  if (is.null(predict)) predict=matrix(mu,nrow=1)
  
  #Computing the upper-triangle matrix of "expectancy of the square"
  v<-apply(#
      apply(predict,1,
            function(pred_i){
              cuhre(ndim=d,ncomp=(d^2+d)/2,
              integrand=function(x){(link.inv(x)%*%t(link.inv(x)))[upper.tri(x%*%t(x),diag=TRUE)]*dmvnorm(x,pred_i,vcov)},
              lower=pred_i-w,upper=pred_i+w,rel.tol=rel.acc,abs.tol= 0.001,
              flags=list(verbose=0))$value}
            ),
      1,mean)
  
  #Creating the VCV matrix
  vcv<-matrix(NA,d,d)
  vcv[upper.tri(vcv,diag=TRUE)]<-v
  vcv[lower.tri(vcv)]<-vcv[upper.tri(vcv)]
  #If necessary, computing the observed multivariate mean
  if (is.null(mvmean.obs))
    mvmean.obs <- QGmvmean(mu,vcov,link.inv,predict=predict,rel.acc=rel.acc,width=width)
  #Computing the VCV matrix using Keonig's formuka
  vcv <- vcv - mvmean.obs%*%t(mvmean.obs)
  
  #Adding the distribution variance if needed (if exp.scale==FALSE)
  if (!exp.scale) {
    vec_vardist = apply(apply(predict,1,function(pred_i){
                  cuhre(ndim=d,ncomp=d,
                  integrand=function(x){var.func(x)*dmvnorm(x,pred_i,vcov)},
                  lower=pred_i-w,upper=pred_i+w,rel.tol=rel.acc,abs.tol=0.0001,
                  flags=list(verbose=0))$value}),1,mean)
    
    vcv <- vcv + diag(vec_vardist)
  }

  #Printing the result
  vcv
}

#Computing the Psi vector
QGmvpsi<-function(mu,vcov,d.link.inv,predict=NULL,rel.acc=0.01,width=10) {
  #Setting the integral width according to vcov (lower mean-w, upper mean+w)
  w<-sqrt(diag(vcov))*width
  #Number of dimensions
  d<-length(w)
  #If predict is not included, then use mu, and 
  if(is.null(predict)) predict=matrix(mu,nrow=1)
  #Computing the mean
  #The double apply is needed to compute the mean for all "predict" values,
  #then average over them
  Psi<-apply(apply(predict,1,function(pred_i){
    cuhre(ndim=d,ncomp=d,
          integrand=function(x){d.link.inv(x)*dmvnorm(x,pred_i,vcov)},
          lower=pred_i-w,upper=pred_i+w,rel.tol=rel.acc,abs.tol=0.0001,
          flags=list(verbose=0))$value}),1,mean)
  #Make Psi a matrix
  Psi<-diag(Psi)
  #print Psi
  Psi
}

##--------------------------------Meta-function for general calculation-----------------------------

QGmvparams<-function(mu,vcv.G,vcv.P,models,predict=NULL,rel.acc=0.01,width=10,n.obs=NULL,theta=NULL,verbose=TRUE)
{
  #Setting the integral width according to vcov (lower mean-w, upper mean+w)
  w<-sqrt(diag(vcv.P))*width
  #Number of dimensions
  d<-length(w)
  #If predict is not included, then use mu, and 
  if(is.null(predict)) predict=matrix(mu,nrow=1)
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
  y_bar<-QGmvmean(mu=mu,vcov=vcv.P,link.inv=inv.links,predict=predict,rel.acc=rel.acc,width=width)
  #Computing the variance-covariance matrix
  if (verbose) print("Computing variance-covariance matrix...")
  vcv.P.obs<-QGvcov(mu=mu,vcov=vcv.P,link.inv=inv.links,var.func=var.funcs,mvmean.obs=y_bar,predict=predict,rel.acc=rel.acc,width=width,exp.scale=FALSE)
  if (verbose) print("Computing Psi...")
  Psi<-QGmvpsi(mu=mu,vcov=vcv.P,d.link.inv=d.inv.links,predict=predict,rel.acc=rel.acc,width=width)
  vcv.G.obs <- Psi %*% vcv.G %*% t(Psi)
  #Return a list of QG parameters on the observed scale
  list(mean.obs=y_bar,vcv.P.obs=vcv.P.obs,vcv.G.obs=vcv.G.obs)
}

##----------------------------------Function to calculate the evolutive prediction-----------------------------

QGmvpred<-function(mu,vcv.G,vcv.P,fit.func,d.fit.func,predict=NULL,rel.acc=0.01,width=10,verbose=TRUE)
{
  #Setting the integral width according to vcov (lower mean-w, upper mean+w)
  w<-sqrt(diag(vcv.P))*width
  #Number of dimensions
  d<-length(w)
  #If predict is not included, then use mu, and 
  if(is.null(predict)) predict=matrix(mu,nrow=1)
  #Calculating the latent mean fitness
  if (verbose) print("Computing mean fitness...")     
  Wbar<-mean(apply(predict,1,function(pred_i){
    cuhre(ndim=d,ncomp=1,                  #Note that ncomp=1 because fitness.func yields a scalar
          integrand=function(x){fit.func(x)*dmvnorm(x,pred_i,vcv.P)},
          lower=pred_i-w,upper=pred_i+w,rel.tol=rel.acc,abs.tol=0.0001,
          flags=list(verbose=0))$value}))
  #Calculating the covariance between latent trait/breeding values and latent fitness
  if (verbose) print("Computing the latent selection and response...")
  #Computing the derivative of fitness
  dW <- apply(apply(predict,1,function(pred_i){
    cuhre(ndim=d,ncomp=d,
          integrand=function(x){d.fit.func(x)*dmvnorm(x,pred_i,vcv.P)},
          lower=pred_i-w,upper=pred_i+w,rel.tol=rel.acc,abs.tol=0.0001,
          flags=list(verbose=0))$value}),1,mean)
  #Computing the selection
  if (dim(predict)[1]>1) sel<-as.vector(((vcv.P+var(predict)) %*% dW)/Wbar) else sel<-as.vector((vcv.P %*% dW)/Wbar)
  #Computing the evolutionary response
  resp<-as.vector((vcv.G %*% dW)/Wbar)
  #Returning the results on the latent scale
  list(mean.lat.fit=Wbar,lat.grad=dW/Wbar,lat.sel=sel,lat.resp=resp)
}