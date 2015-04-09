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
QGmean.obs.mv<-function(mu,vcov,link.inv,predict=NULL,rel.acc=0.01,width=10) {
  #Setting the integral width according to vcov (lower mean-w, upper mean+w)
  w<-sqrt(diag(vcov))*width
  #Number of dimensions
  d<-length(w)
  #If predict is not included, then use mu, and 
  if(is.null(predict)) predict=matrix(mu,nrow=1)
  #Computing the mean
  apply(apply(predict,1,function(pred_i){
      cuhre(ndim=d,ncomp=d,
            integrand=function(x){link.inv(x)*dmvnorm(x,pred_i,vcov)},
            lower=pred_i-w,upper=pred_i+w,rel.tol=rel.acc,abs.tol=0.0001,
            flags=list(verbose=0))$value}),1,mean)
}