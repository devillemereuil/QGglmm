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

##---------------------------------General functions----------------------------------------

#Calculating the observed/expected scale mean
QGmean.obs<-function(mu,var,link.inv,width=35,predict=NULL) {
  #If no fixed effects were included in the model
  if (is.null(predict)) predict=0
  mean(sapply(predict,function(pred_i){integrate(f=function(x){link.inv(x)*dnorm(x,mu+pred_i,sqrt(var))},lower=mu+pred_i-width*sqrt(var),upper=mu+pred_i+width*sqrt(var))$value}))
}

#Calculating the expected scale variance
QGvar.exp<-function(mu,var,link.inv,obs.mean=NULL,width=35,predict=NULL) {
  #If not provided, compute the obsereved mean
  if (is.null(obs.mean)){
    obs.mean=QGmean.obs(mu=mu,var=var,link.inv=link.inv,width=width,predict=predict)
  }
  if (is.null(predict)) predict=0
  #Note: Using Koenig's formula
  mean(sapply(predict,function(pred_i){integrate(f=function(x){(link.inv(x)**2)*dnorm(x,mu+pred_i,sqrt(var))},lower=mu+pred_i-width*sqrt(var),upper=mu+pred_i+width*sqrt(var))$value}))-(obs.mean**2)
}

#Calculating the "distribution" variance
QGvar.dist<-function(mu,var,var.func,width=35,predict=NULL) {
  #If no fixed effects were included in the model
  if (is.null(predict)) predict=0
  mean(sapply(predict,function(pred_i){integrate(f=function(x){var.func(x)*dnorm(x,mu+pred_i,sqrt(var))},lower=mu+pred_i-width*sqrt(var),upper=mu+pred_i+width*sqrt(var))$value}))
}

#Calculating "psi" for the observed additive genetic variance computation
QGpsi<-function(mu,var,d.link.inv,width=35,predict=NULL) {
    #If no fixed effects were included in the model
  if (is.null(predict)) predict=0
    psi<-integrate(f=function(x){d.link.inv(x)*dnorm(x,mu,sqrt(var))},lower=mu-width*sqrt(var),upper=mu+width*sqrt(var))$value
    mean(sapply(predict,function(pred_i){integrate(f=function(x){d.link.inv(x)*dnorm(x,mu+pred_i,sqrt(var))},lower=mu+pred_i-width*sqrt(var),upper=mu+pred_i+width*sqrt(var))$value}))
}

##-------------------------------"Dictionary" functions---------------------------------------------

#Function creating the needed functions according to the link name
QGlink.funcs<-function(name,n.obs=NULL,theta=NULL) {
  if (name=="binom1.probit") {
    inv.link=function(x){pnorm(x)}
    var.func=function(x){pnorm(x)*(1-pnorm(x))}
    d.inv.link=function(x){dnorm(x)}
  } else  if (name=="binomN.probit") {
    if (is.null(n.obs)) {stop("binomN model used, but no observation number (n.obs) defined.")}
    inv.link=function(x){n.obs*pnorm(x)}
    var.func=function(x){n.obs*pnorm(x)*(1-pnorm(x))}
    d.inv.link=function(x){n.obs*dnorm(x)}
  } else if (name=="binom.logit") {
    inv.link=function(x){plogis(x)}
    var.func=function(x){plogis(x)*(1-plogis(x))}
    d.inv.link=function(x){dlogis(x)}
  } else if (name=="binomN.logit") {
    if (is.null(n.obs)) {stop("binomN model used, but no observation number (n.obs) defined.")}
    inv.link=function(x){n.obs*plogis(x)}
    var.func=function(x){n.obs*plogis(x)*(1-plogis(x))}
    d.inv.link=function(x){n.obs*dlogis(x)}
  } else if (name=="threshold") {
      ##TODO
      stop("Not implemented yet")
  } else if (name=="Poisson.log") {
    inv.link=function(x){exp(x)}
    var.func=function(x){exp(x)}
    d.inv.link=function(x){exp(x)}
  } else if (name=="Poisson.sqrt") {
    inv.link=function(x){x**2}
    var.func=function(x){x**2}
    d.inv.link=function(x){2*x}
  } else if (name=="negbin.log") {
    if (is.null(theta)) {stop("negbin model used, but theta not defined.")}
    inv.link=function(x){exp(x)}
    var.func=function(x){exp(x)+(exp(2*x)/theta)}
    d.inv.link=function(x){exp(x)}
  } else if (name=="negbin.sqrt") {
    if (is.null(theta)) {stop("negbin model used, but theta not defined.")}
    inv.link=function(x){x**2}
    var.func=function(x){(x**2)+((x**4)/theta)}
    d.inv.link=function(x){2*x}
  } else {stop("Invalid model name. Use a valid model name or enter a custom model specification.")}
  list(inv.link=inv.link,var.func=var.func,d.inv.link=d.inv.link)
}

##-----------------------------Special functions for known analytical solutions--------------

qg.binom1.probit=function(mu,var.a,var.p,predict=NULL) {
  if (is.null(predict)) predict=0;
  #Observed mean
  p=mean(1-pnorm(0,mu+predict,sqrt(var.p+1)))
  #Observed variance
  var_obs=p*(1-p)
  #Psi
  Psi=mean(dnorm(0,(mu+predict),sqrt(var.p+1)))
  data.frame(mean.obs=p,var.obs=var_obs,var.a.obs=(Psi**2)*var.a,h2.obs=((Psi**2)*var.a)/var_obs)
}

qg.binomN.probit=function(mu,var.a,var.p,n.obs,predict=NULL,width=35) {
  if (is.null(predict)) predict=0;
  #Observed mean
  p=n.obs*mean(1-pnorm(0,mu+predict,sqrt(var.p+1)))
  #Observed variance
  prob.sq.int=mean(sapply(predict,function(pred_i){integrate(f=function(x){(pnorm(x)**2)*dnorm(x,mu+pred_i,sqrt(var.p))},lower=mu+pred_i-width*sqrt(var.p),upper=mu+pred_i+width*sqrt(var.p))$value}))
  var_obs=((n.obs**2)-n.obs)*prob.sq.int - p**2 +p
  #Psi
  Psi=n.obs*mean(dnorm(0,(mu+predict),sqrt(var.p+1)))
  data.frame(mean.obs=p,var.obs=var_obs,var.a.obs=(Psi**2)*var.a,h2.obs=((Psi**2)*var.a)/var_obs)
}

qg.Poisson.log=function(mu,var.a,var.p,predict=NULL) {
  if (is.null(predict)) predict=0
  #Observed mean
  lambda=mean(exp(mu+predict+(var.p/2)))
  #Mean of lambda square, needed for the following
  lambda_sq=mean(exp(2*(mu+predict+var.p/2)))
  #Observed variance
  var_obs=lambda_sq*exp(var.p)-lambda**2+lambda
  data.frame(mean.obs=lambda,var.obs=var_obs,var.a.obs=(lambda**2)*var.a,h2.obs=((lambda**2)*var.a)/var_obs)
}

qg.Poisson.sqrt=function(mu,var.a,var.p,predict=NULL) {
  if (is.null(predict)) predict=0
  #Observed mean
  lambda=mean((mu+predict)**2+var.p)
  #Observed variance
  var_obs=mean((mu+predict)**4 + 6*var.p*((mu+predict)**2)+3*(var.p**2))-lambda**2+lambda
  #Psi
  Psi=mean(2*(mu+predict))
  data.frame(mean.obs=lambda,var.obs=var_obs,var.a.obs=(Psi**2)*var.a,h2.obs=((Psi**2)*var.a)/var_obs)
}

qg.negbin.log=function(mu,var.a,var.p,theta,predict=NULL) {
  if (is.null(predict)) predict=0
  #Observed mean
  lambda=mean(exp(mu+predict+(var.p/2)))
  #Mean of lambda square, needed for the following
  lambda_sq=mean(exp(2*(mu+predict+var.p/2)))
  #Observed variance
  var_obs=lambda_sq*exp(var.p)-lambda**2+lambda+mean(exp(2*(mu+predict+var.p)))/theta
  data.frame(mean.obs=lambda,var.obs=var_obs,var.a.obs=(lambda**2)*var.a,h2.obs=((lambda**2)*var.a)/var_obs)
}

qg.negbin.sqrt=function(mu,var.a,var.p,theta,predict=NULL) {
  if (is.null(predict)) predict=0
  #Observed mean
  lambda=mean((mu+predict)**2+var.p)
  #Observed variance
  var_obs=mean((mu+predict)**4 + 6*var.p*((mu+predict)**2)+3*(var.p**2))-lambda**2+lambda+(mean((mu+predict)**4 + 6*var.p*((mu+predict)**2)+3*(var.p**2))/theta)
  #Psi
  Psi=mean(2*(mu+predict))
  data.frame(mean.obs=lambda,var.obs=var_obs,var.a.obs=(Psi**2)*var.a,h2.obs=((Psi**2)*var.a)/var_obs)
}

##--------------------------------Meta-function for general calculation-----------------------------

QGparams<-function(mu,var.a,var.p,model="",width=35,predict=NULL,closed.form=TRUE,custom.model=NULL,n.obs=NULL,theta=NULL,verbose=TRUE) {

  ##Using analytical solutions if possible (and asked for, see closed.form arg)
  if (model=="binom1.probit"&closed.form) {						#Binary.probit model
      if (verbose) print("Using the closed forms for a Binomial1-probit model.")
      qg.binom1.probit(mu=mu,var.a=var.a,var.p=var.p,predict=predict)
  } else if (model=="binomN.probit"&closed.form) {					#Binomial-not-binary model
      if (is.null(n.obs)) {stop("binomN.probit model used, but no observation number (n.obs) defined.")}
      if (verbose) print("Using the closed forms for a BinomialN-probit model.")
      qg.binomN.probit(mu=mu,var.a=var.a,var.p=var.p,predict=predict,n.obs=n.obs,width=width)
  } else if (model=="Poisson.log"&closed.form){						#Poisson-log model
      if(verbose) print("Using the closed forms for a Poisson-log model.")
      qg.Poisson.log(mu=mu,var.a=var.a,var.p=var.p,predict=predict)
  } else if (model=="Poisson.sqrt"&closed.form){						#Poisson-sqrt model
      if(verbose) print("Using the closed forms for a Poisson-sqrt model.")
      qg.Poisson.sqrt(mu=mu,var.a=var.a,var.p=var.p,predict=predict)
  } else if (model=="negbin.log"&closed.form){						#NegBin-log model
      if (is.null(theta)) {stop("negbin model used, but theta not defined.")}
      if(verbose) print("Using the closed forms for a NegativeBinomial-log model.")
      qg.negbin.log(mu=mu,var.a=var.a,var.p=var.p,predict=predict,theta=theta)
  } else if (model=="negbin.sqrt"&closed.form){						#NegBin-sqrt model
      if (is.null(theta)) {stop("negbin model used, but theta not defined.")}
      if(verbose) print("Using the closed forms for a NegativeBinomial-sqrt model.")
      qg.negbin.sqrt(mu=mu,var.a=var.a,var.p=var.p,predict=predict,theta=theta)
  } else {
  
  ##Else, use the general integral equations
  #Use a custom model if defined, otherwise look into the "Dictionary"
      if (is.null(custom.model)) {
        if (model==""){stop("The function requires either model or custom.model.")} else {
          funcs=QGlink.funcs(model,n.obs=n.obs,theta=theta)
        }} else {funcs=custom.model}
  #Observed mean computation
      if (verbose) print("Computing observed mean...")
      y_bar=QGmean.obs(mu=mu,var=var.p,link.inv=funcs$inv.link,width=width,predict=predict)
  #Variances computation
      if (verbose) print("Computing variances...")
      var_exp=QGvar.exp(mu=mu,var=var.p,link.inv=funcs$inv.link,obs.mean=y_bar,width=width,predict=predict)
      var_dist=QGvar.dist(mu=mu,var=var.p,var.func=funcs$var.func,width=width,predict=predict)
      var_obs=var_exp+var_dist
  #Psi computation (for the observed additive genetic variance)
      if (verbose) print("Computing Psi...")
      Psi=QGpsi(mu=mu,var=var.p,d.link.inv=funcs$d.inv.link,width=width,predict=predict)
  #Return a data.frame with the calculated components
      data.frame(mean.obs=y_bar,var.obs=var_obs,var.a.obs=(Psi**2)*var.a,h2.obs=((Psi**2)*var.a)/var_obs)
  }
}

##----------------------------------Function to calculate the evolutive prediction-----------------------------

QGpred<-function(mu,var.a,var.p,model,fitness.func,width=35,predict=NULL,custom.model=NULL,n.obs=NULL,theta=NULL,verbose=TRUE) {
  if (is.null(predict)) predict=0
  #Calculating the latent mean fitness
  print("Computing mean fitness...")
  Wbar=mean(sapply(predict,function(pred_i){integrate(f=function(x){fitness.func(x)*dnorm(x,mu+pred_i,sqrt(var.p))},lower=mu+pred_i-width*sqrt(var.p),upper=mu+pred_i+width*sqrt(var.p))$value}))
  #Calculating the covariance between latent trait and latent fitness
  print("Computing the latent selection...")
  sel=(mean(sapply(predict,function(pred_i){integrate(f=function(x){x*fitness.func(x)*dnorm(x,mu+pred_i,sqrt(var.p))},lower=mu+pred_i-width*sqrt(var.p),upper=mu+pred_i+width*sqrt(var.p))$value}))-(mean(mu+predict)*Wbar))/Wbar
  print(sel)
  #Calculating the covariance between the breeding values and the latent fitness
  print("Computing the latent response... (this might take a while if predict is large)")
  #Calculating the latent mean fitness conditional to a
  exp_W_a<-function(vec){sapply(vec,function(a){mean(sapply(predict,function(pred_i){integrate(f=function(x){fitness.func(x)*dnorm(x,mu+a+pred_i,sqrt(var.p-var.a))},lower=mu+a+pred_i-width*sqrt(var.p-var.a),upper=mu+a+pred_i+width*sqrt(var.p-var.a))$value}))})}
  resp=(integrate(f=function(a){a*exp_W_a(a)*dnorm(a,0,sqrt(var.a))},lower=-width*sqrt(var.a),upper=width*sqrt(var.a))$value)/Wbar
  #Returning the results on the latent scale
  data.frame(mean.lat.fit=Wbar,pred.lat.sel=sel,pred.lat.resp=resp)
}