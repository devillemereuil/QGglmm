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
QGmvmean.obs<-function(mu,vcov,link.inv,predict=NULL,rel.acc=0.01,width=10) {
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
QGvcov.obs<-function(mu,vcov,link.inv,var.func,mvmean.obs=NULL,predict=NULL,rel.acc=0.01,width=10,exp.scale=FALSE) {
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
    mvmean.obs <- QGmvmean.obs(mu,vcov,link.inv,predict=predict,rel.acc=rel.acc,width=width)
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
QGpsi.mv<-function(mu,vcov,d.link.inv,predict=NULL,rel.acc=0.01,width=10) {
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
          integrand=function(x){d.link.inv(x)*dmvnorm(x,pred_i,vcov)},
          lower=pred_i-w,upper=pred_i+w,rel.tol=rel.acc,abs.tol=0.0001,
          flags=list(verbose=0))$value}),1,mean)
}