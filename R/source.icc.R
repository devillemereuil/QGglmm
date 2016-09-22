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

QGicc <- function(mu, var.comp, var.p, link.inv, predict=NULL,width=10){
  if(length(mu)>1 | length(var) >1) stop("The parameters mu and var must be of length 1, please check your input.")
  if (is.null(predict)) { if(is.null(mu)) {stop("Please provide either mu or predict.")} else {predict=mu;}}
  
  error("Nothing for now!")
}


QGmvicc<-function(mu,vcov.comp,vcov.p,link.inv,mvmean.obs=NULL,predict=NULL,rel.acc=0.01,width=10) {
  #Setting the integral width according to vcov (lower mean-w, upper mean+w)
  w<-sqrt(diag(vcov.p))*width
  #Number of dimensions
  d<-length(w)
  #If no fixed effects were included in the model
  if (is.null(predict)) predict=matrix(mu,nrow=1)
  
  #Function giving the conditional expectancy
  ##TODO: make this function predict-aware (i.e. replace mu by pred_i and average over it)
  cond_exp <- function(a) {
    cuhre(ndim=d,ncomp=d,
              integrand=function(x){link.inv(x)*dmvnorm(x,mu+a,vcov.p-vcov.comp)},
              lower=mu-w,upper=mu+w,rel.tol=rel.acc,abs.tol= 0.001,
              flags=list(verbose=0))$value
  }
  
  #Computing the upper-triangle matrix of "expectancy of the square"
  v<-apply(#
      apply(predict,1,
            function(pred_i){
              cuhre(ndim=d,ncomp=(d^2+d)/2,
              integrand=function(x){
		  (cond_exp(x)%*%t(cond_exp(x)))[upper.tri(x%*%t(x),diag=TRUE)]*dmvnorm(x,rep(0,d),vcov.comp)
	      },
	      lower=-w,upper=w,rel.tol=rel.acc,abs.tol= 0.001,
	      flags=list(verbose=0))$value
	    }
      ),
     1,mean)
  
  #Creating the VCV matrix
  vcv<-matrix(NA,d,d)
  vcv[upper.tri(vcv,diag=TRUE)]<-v
  vcv[lower.tri(vcv)]<-vcv[upper.tri(vcv)]
  #If necessary, computing the observed multivariate mean
  if (is.null(mvmean.obs))
    mvmean.obs <- QGmvmean(mu,vcov.p,link.inv,predict=predict,rel.acc=rel.acc,width=width)
  #Computing the VCV matrix using Keonig's formuka
  vcv <- vcv - mvmean.obs%*%t(mvmean.obs)
  
  #Printing the result
  vcv
} 
