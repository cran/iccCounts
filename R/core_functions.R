

#' Estimates GLMM. Concordance setting.
#'
#' Estimates the GLMM model for concordance setting
#' @keywords internal
#' @param data A data frame containing at least three columns: outcome (named as y), subject identifier (named as id) and method identifier (named as met).
#' @param fam Character string. The within-subjects pdf to use. Valid options are: "poisson" (default) for Poisson pdf; "nbinom1" for Negative Binomial pdf with additive extradispersion; "nbinom2" for Negative Binomial pdf with proportional extradispersion; "zip" for zero-inflated Poisson pdf; "zinb1" for zero-inflated nbinom1 pdf; "zinb2" for zero-inflated nbinom2 pdf;
#' @return An object of class glmmTMB with the model estimates.

fit_model_con<-function(data,fam=c("poisson","nbinom2","nbinom1","zip","zinb1","zinb2")){

  fam<-match.arg(fam)
  if (length(data[,"met"])==0) stop("No methods variable")

  if (fam=="poisson") model<- glmmTMB(y~met+(1|id),data=data,family=poisson)
  if (fam=="nbinom2") model<- glmmTMB(y~met+(1|id),data=data,family=nbinom2())
  if (fam=="nbinom1") model<- glmmTMB(y~met+(1|id),data=data,family=nbinom1())
  if (fam=="zip") model<- glmmTMB(y~met+(1|id),data=data,ziformula = ~1,family=poisson())
  if (fam=="zinb1") model<- glmmTMB(y~met+(1|id),data=data,ziformula = ~1,family=nbinom1())
  if (fam=="zinb2") model<- glmmTMB(y~met+(1|id),data=data,ziformula = ~1,family=nbinom2())

  return(model)


}

#' Estimates GLMM. Repeatability setting.
#'
#' Estimates the GLMM model for concordance setting
#' @keywords internal
#' @param data A data frame containing at least two columns: outcome (named as y) and subject identifier (named as id).
#' @param fam Character string. The within-subjects pdf to use. Valid options are: "poisson" (default) for Poisson pdf; "nbinom1" for Negative Binomial pdf with additive extradispersion; "nbinom2" for Negative Binomial pdf with proportional extradispersion; "zip" for zero-inflated Poisson pdf; "zinb1" for zero-inflated nbinom1 pdf; "zinb2" for zero-inflated nbinom2 pdf;
#' @return An object of class glmmTMB with the model estimates.

fit_model_rep<-function(data,fam=c("poisson","nbinom2","nbinom1","zip","zinb1","zinb2")){

  fam<-match.arg(fam,choices=c("poisson","nbinom2","nbinom1","zip","zinb1","zinb2"))

  if (fam=="poisson") model<- glmmTMB(y~(1|id),data=data,family=poisson)
  if (fam=="nbinom2") model<- glmmTMB(y~(1|id),data=data,family=nbinom2())
  if (fam=="nbinom1") model<- glmmTMB(y~(1|id),data=data,family=nbinom1())
  if (fam=="zip") model<- glmmTMB(y~(1|id),data=data,ziformula = ~1,family=poisson())
  if (fam=="zinb1") model<- glmmTMB(y~(1|id),data=data,ziformula = ~1,family=nbinom1())
  if (fam=="zinb2") model<- glmmTMB(y~(1|id),data=data,ziformula = ~1,family=nbinom2())

  return(model)


}

#' ICC concordance
#'
#' Estimates the ICC for concordance setting for poisson, nbinom1 and nbinom2 families.
#' @keywords internal
#' @param data A data frame containing at least three columns: outcome (named as y), subject identifier (named as id) and method identifier (named as met).
#' @param model A glmmTMB obejct.
#' @return A data frame with ICC estimate, its standard error and variance components.

r_est_con<-function(data,model){

  # Model estimates
  pars<-model$fit$par

  # between-observers variability
  k<-length(unique(data$met))
  if (k<2) stop("Less than two methods")
  b<-pars[1:k]
  sb<-sb_f(k,b,varB(model,k))

  fam<-model$modelInfo$family$family

  if (fam=="poisson"){
    r_est<-r_Pois2(pars[1],pars[k+1],sb)
    CovB<-model$sdr$cov.fixed[c(1,k+1),c(1,k+1)]
    se.r<-se_r_Pois2(pars[1],pars[k+1],sb,CovB,var_sb(k,b,varB(model,k)))
    sa<-sa_f(pars[k+1])
    mu<-mu_f2(pars[1],pars[k+1],sb)
    out<-data.frame(r_est,se.r,mu,sa,sb)
    names(out)<-c("ICC","SE ICC","mu","BSVar","BMVar")

  }

  if (fam=="nbinom2"){
    r_est<-r_NB2_2(pars[1],pars[k+2],pars[k+1],sb)
    CovB<-model$sdr$cov.fixed[c(1,k+2,k+1),c(1,k+2,k+1)]
    se.r<-se_r_NB2_2(pars[1],pars[k+2],pars[k+1],sb,CovB,var_sb(k,b,varB(model,k)))
    sa<-sa_f(pars[k+2])
    mu<-mu_f2(pars[1],pars[k+2],sb)
    r<-1/k_f(pars[k+1])
    out<-data.frame(r_est,se.r,mu,sa,sb,r)
    names(out)<-c("ICC","SE ICC","mu","BSVar","BMVar","r")

  }

  if (fam=="nbinom1"){
    r_est<-r_NB1_2(pars[1],pars[k+2],pars[k+1],sb)
    CovB<-model$sdr$cov.fixed[c(1,k+2,k+1),c(1,k+2,k+1)]
    se.r<-se_r_NB1_2(pars[1],pars[k+2],pars[k+1],sb,CovB,var_sb(k,b,varB(model,k)))
    sa<-sa_f(pars[k+2])
    mu<-mu_f2(pars[1],pars[k+2],sb)
    r<-k_f(pars[k+1])
    out<-data.frame(r_est,se.r,mu,sa,sb,r)
    names(out)<-c("ICC","SE ICC","mu","BSVar","BMVar","r")

  }

  return(out)

}

#' ICC concordance for zero-inflated models
#'
#' Estimates the ICC for concordance setting for zero-inflated models (poisson, nbinom1 and nbinom2 families).
#' @keywords internal
#' @param data A data frame containing at least three columns: outcome (named as y), subject identifier (named as id) and method identifier (named as met).
#' @param model A glmmTMB obejct.
#' @return A data frame with ICC estimate, its standard error and variance components.

r_est_con_zi<-function(data,model){

  # Model estimates
  pars<-model$fit$par

  k<-length(unique(data$met))
  if (k<2) stop("Less than two methods")
  b<-pars[1:k]
  sb<-sb_f(k,b,varB(model,k))
  fam<-model$modelInfo$family$family

  if (fam=="poisson"){
    r_est<-r_ZIP2(pars[1],pars[k+2],pars[k+1],sb)
    CovB<-model$sdr$cov.fixed[c(1,k+2,k+1),c(1,k+2,k+1)]
    se.r<-se_r_ZIP2(pars[1],pars[k+2],pars[k+1],sb,CovB,var_sb(k,b,varB(model,k)))
    mu<-mu_f2(pars[1],pars[k+2],sb)
    sa<-sa_f(pars[k+2])
    pi<-p_f(pars[k+1])
    out<-data.frame(r_est,se.r,mu,sa,sb,pi)
    names(out)<-c("ICC","SE ICC","mu","BSVar","BMVar","Pi")
  }

  if (fam=="nbinom2"){
    r_est<-r_ZINB2_2(pars[1],pars[k+3],pars[k+2],pars[k+1],sb)
    CovB<-model$sdr$cov.fixed[c(1,k+3,k+2,k+1),c(1,k+3,k+2,k+1)]
    se.r<-se_r_ZINB2_2(pars[1],pars[k+3],pars[k+2],pars[k+1],sb,CovB,var_sb(k,b,varB(model,k)))
    mu<-mu_f2(pars[1],pars[k+3],sb)
    sa<-sa_f(pars[k+3])
    r<-1/k_f(pars[k+2])
    pi<-p_f(pars[k+1])
    varcomp<-c(mu,sa,sb,r,pi)
    out<-data.frame(r_est,se.r,mu,sa,sb,r,pi)
    names(out)<-c("ICC","SE ICC","mu","BSVar","BMVar","r","pi")
  }

  if (fam=="nbinom1"){
    r_est<-r_ZINB1_2(pars[1],pars[k+3],pars[k+2],pars[k+1],sb)
    CovB<-model$sdr$cov.fixed[c(1,k+3,k+2,k+1),c(1,k+3,k+2,k+1)]
    se.r<-se_r_ZINB1_2(pars[1],pars[k+3],pars[k+2],pars[k+1],sb,CovB,var_sb(k,b,varB(model,k)))
    mu<-mu_f2(pars[1],pars[k+3],sb)
    sa<-sa_f(pars[k+3])
    r<-k_f(pars[k+2])
    pi<-p_f(pars[k+1])
    varcomp<-c(mu,sa,sb,r,pi)
    out<-data.frame(r_est,se.r,mu,sa,sb,r,pi)
    names(out)<-c("ICC","SE ICC","mu","BSVar","BMVar","r","pi")

  }



  return(out)

}


#' ICC repeatability
#'
#' Estimates the ICC for repeatability setting for poisson, nbinom1 and nbinom2 families.
#' @keywords internal
#' @param data A data frame containing at least two columns: outcome (named as y) and subject identifier (named as id).
#' @param model A glmmTMB obejct.
#' @return A data frame with ICC estimate, its standard error and variance components.


r_est_rep<-function(data,model){

  # Model estimates
  pars<-model$fit$par

  fam<-model$modelInfo$family$family

  if (fam=="poisson"){
    r_est<-r_Pois(pars[1],pars[2])
    CovB<-model$sdr$cov.fixed
    se.r<-se_r_Pois(pars[1],pars[2],CovB)
    mu<-mu_f(pars[1],pars[2])
    sa<-sa_f(pars[2])
    out<-data.frame(r_est,se.r,mu,sa)
    names(out)<-c("ICC","SE ICC","mu","BSVar")
  }

  if (fam=="nbinom2"){
    r_est<-r_NB2(pars[1],pars[3],pars[2])
    CovB<-model$sdr$cov.fixed[c(1,3,2),c(1,3,2)]
    se.r<-se_r_NB2(pars[1],pars[3],pars[2],CovB)
    mu<-mu_f(pars[1],pars[3])
    sa<-sa_f(pars[3])
    r<-1/k_f(pars[2])
    varcomp<-c(mu,sa,r)
    out<-data.frame(r_est,se.r,mu,sa,r)
    names(out)<-c("ICC","SE ICC","mu","BSVar","r")

  }

  if (fam=="nbinom1"){
    r_est<-r_NB1(pars[1],pars[3],pars[2])
    CovB<-model$sdr$cov.fixed[c(1,3,2),c(1,3,2)]
    se.r<-se_r_NB1(pars[1],pars[3],pars[2],CovB)
    mu<-mu_f(pars[1],pars[3])
    sa<-sa_f(pars[3])
    r<-k_f(pars[2])
    varcomp<-c(mu,sa,r)
    out<-data.frame(r_est,se.r,mu,sa,r)
    names(out)<-c("ICC","SE ICC","mu","BSVar","r")

  }



  return(out)

}

#' ICC repeatability for zero-inflated models
#'
#' Estimates the ICC for repeatability setting for zero-inflated models (poisson, nbinom1 and nbinom2 families).
#' @keywords internal
#' @param data A data frame containing at least two columns: outcome (named as y) and subject identifier (named as id).
#' @param model A glmmTMB obejct.
#' @return A data frame with ICC estimate, its standard error and variance components.


r_est_rep_zi<-function(data,model){

  # Model estimates
  pars<-model$fit$par

  fam<-model$modelInfo$family$family

  if (fam=="poisson"){
    r_est<-r_ZIP(pars[1],pars[3],pars[2])
    CovB<-model$sdr$cov.fixed[c(1,3,2),c(1,3,2)]
    se.r<-se_r_ZIP(pars[1],pars[3],pars[2],CovB)
    mu<-mu_f(pars[1],pars[3])
    sa<-sa_f(pars[3])
    pi<-p_f(pars[2])
    out<-data.frame(r_est,se.r,mu,sa,pi)
    names(out)<-c("ICC","SE ICC","mu","BSVar","Pi")
  }

  if (fam=="nbinom2"){
    r_est<-r_ZINB2(pars[1],pars[4],pars[3],pars[2])
    CovB<-model$sdr$cov.fixed[c(1,4,3,2),c(1,4,3,2)]
    se.r<-se_r_ZINB2(pars[1],pars[4],pars[3],pars[2],CovB)
    mu<-mu_f(pars[1],pars[4])
    sa<-sa_f(pars[4])
    pi<-p_f(pars[2])
    r<-1/k_f(pars[3])
    out<-data.frame(r_est,se.r,mu,sa,r,pi)
    names(out)<-c("ICC","SE ICC","mu","BSVar","r","Pi")
  }

  if (fam=="nbinom1"){
    r_est<-r_ZINB1(pars[1],pars[4],pars[3],pars[2])
    CovB<-model$sdr$cov.fixed[c(1,4,3,2),c(1,4,3,2)]
    se.r<-se_r_ZINB1(pars[1],pars[4],pars[3],pars[2],CovB)
    mu<-mu_f(pars[1],pars[4])
    sa<-sa_f(pars[4])
    pi<-p_f(pars[2])
    r<-k_f(pars[3])
    out<-data.frame(r_est,se.r,mu,sa,r,pi)
    names(out)<-c("ICC","SE ICC","mu","BSVar","r","Pi")

  }



  return(out)

}

#'  ICC Poisson model. Repeatability setting
#'
#'  ICC for Poisson model. Repeatability setting
#' @keywords internal
#' @param b1 Intercept of GLMM
#' @param b2 Parameter related to between-subjects variance in the GLMM model
#' @return Scalar

r_Pois<-function(b1,b2){
  mu_f(b1,b2)*(exp(sa_f(b2))-1)/(mu_f(b1,b2)*(exp(sa_f(b2))-1)+1)
}


#'  ICC Poisson model. Concordance setting
#'
#'  ICC for Poisson model. Concordance setting
#' @keywords internal
#' @param b1 Intercept of GLMM
#' @param b2 Parameter related to between-subjects variance in the GLMM model
#' @param sb Between-methods variability
#' @return Scalar

r_Pois2<-function(b1,b2,sb){
  mu_f2(b1,b2,sb)*(exp(sa_f(b2))-1)/(mu_f2(b1,b2,sb)*(exp(sa_f(b2)+sb)-1)+1)
}

#'  ICC standard error for Poisson model. Repeatability setting
#'
#'  ICC standard error for Poisson model. Repeatability setting
#' @keywords internal
#' @param b1 Intercept of GLMM
#' @param b2 Parameter related to between-subjects variance in the GLMM model
#' @param CovB Variance covariance matrix of b1 and b2
#'
#' @return Scalar

se_r_Pois<-function(b1,b2,CovB){

  d1<-Deriv(r_Pois,"b1")
  d2<-Deriv(r_Pois,"b2")

  dev.theta<-matrix(
    c(d1(b1,b2) ,
      d2(b1,b2))
    ,nrow=2)

  sqrt(t(dev.theta)%*%CovB%*%dev.theta)

}

#'  ICC standard error for Poisson model. Concordance setting
#'
#'  ICC standard error for Poisson model. Concordance setting
#' @keywords internal
#' @param b1 Intercept of GLMM
#' @param b2 Parameter related to between-subjects variance in the GLMM model
#' @param sb Between-methods variability
#' @param CovB Variance covariance matrix of b1 and b2
#' @param var.sb Standard error of sb (in variance scale).
#' @return Scalar

se_r_Pois2<-function(b1,b2,sb,CovB,var.sb){

  d1<-Deriv(r_Pois2,"b1")
  d2<-Deriv(r_Pois2,"b2")
  d3<-Deriv(r_Pois2,"sb")

  dev.theta<-matrix(
    c(d1(b1,b2,sb) ,
      d2(b1,b2,sb),
      d3(b1,b2,sb))
    ,nrow=3)


  S<-array(0,c(3,3))
  S[1:2,1:2]<-CovB
  S[3,3]<-var.sb

  sqrt(t(dev.theta)%*%S%*%dev.theta)

}

#'  ICC for Negative Binomial with proportional extradispersion. Repeatability setting
#'
#'  ICC for Negative Binomial with proportional extradispersion. Repeatability setting
#' @keywords internal
#' @param b1 Intercept of GLMM
#' @param b2 Parameter related to between-subjects variance in the GLMM model
#' @param b3 Parameter related to extra-dispersion in the GLMM model
#'
#' @return Scalar


r_NB2<-function(b1,b2,b3){
  K<-(1+k_f(b3))/(k_f(b3))
  mu_f(b1,b2)*(exp(sa_f(b2))-1)/(mu_f(b1,b2)*(K*exp(sa_f(b2))-1)+1)
}

#'  ICC standard error for Negative Binomial with proportional extradispersion. Repeatability setting
#'
#'  ICC standard error for Negative Binomial with proportional extradispersion. Repeatability setting
#' @keywords internal
#' @param b1 Intercept of GLMM
#' @param b2 Parameter related to between-subjects variance in the GLMM model
#' @param b3 Parameter related to extra-dispersion in the GLMM model
#' @param CovB Variance covariance matrix of b1, b2 and b3
#'
#' @return Scalar

se_r_NB2<-function(b1,b2,b3,CovB){

  d1<-Deriv(r_NB2,"b1")
  d2<-Deriv(r_NB2,"b2")
  d3<-Deriv(r_NB2,"b3")

  dev.theta<-matrix(
    c(d1(b1,b2,b3) ,
      d2(b1,b2,b3),
      d3(b1,b2,b3))
    ,nrow=3)


  S<-array(0,c(3,3))
  S[1:3,1:3]<-CovB

  sqrt(t(dev.theta)%*%S%*%dev.theta)

}

#'  ICC for Negative Binomial with proportional extradispersion. Concordance setting
#'
#'  ICC for Negative Binomial with proportional extradispersion. Concordance setting
#' @keywords internal
#' @param b1 Intercept of GLMM
#' @param b2 Parameter related to between-subjects variance in the GLMM model
#' @param b3 Parameter related to extra-dispersion in the GLMM model
#' @param sb Between-methods variability
#' @return Scalar

r_NB2_2<-function(b1,b2,b3,sb){
  K<-(1+k_f(b3))/(k_f(b3))
  mu_f2(b1,b2,sb)*(exp(sa_f(b2))-1)/(mu_f2(b1,b2,sb)*(K*exp(sa_f(b2)+sb)-1)+1)
}

#'  ICC standard error for Negative Binomial with proportional extradispersion. Repeatability setting
#'
#'  ICC standard error for Negative Binomial with proportional extradispersion. Repeatability setting
#' @keywords internal
#' @param b1 Intercept of GLMM
#' @param b2 Parameter related to between-subjects variance in the GLMM model
#' @param b3 Parameter related to extra-dispersion in the GLMM model
#' @param sb Between-methods variability
#' @param CovB Variance covariance matrix of b1, b2 and b3
#' @param var.sb Standard error of sb (in variance scale).
#' @return Scalar
se_r_NB2_2<-function(b1,b2,b3,sb,CovB,var.sb){

  d1<-Deriv(r_NB2_2,"b1")
  d2<-Deriv(r_NB2_2,"b2")
  d3<-Deriv(r_NB2_2,"b3")
  d4<-Deriv(r_NB2_2,"sb")

  dev.theta<-matrix(
    c(d1(b1,b2,b3,sb) ,
      d2(b1,b2,b3,sb),
      d3(b1,b2,b3,sb),
      d4(b1,b2,b3,sb))
    ,nrow=4)


  S<-array(0,c(4,4))
  S[1:3,1:3]<-CovB
  S[4,4]<-var.sb

  sqrt(t(dev.theta)%*%S%*%dev.theta)

}

#'  ICC for Negative Binomial with additive extradispersion. Repeatability setting
#'
#'  ICC for Negative Binomial with additive extradispersion. Repeatability setting
#' @keywords internal
#' @param b1 Intercept of GLMM
#' @param b2 Parameter related to between-subjects variance in the GLMM model
#' @param b3 Parameter related to extra-dispersion in the GLMM model
#' @return Scalar

r_NB1<-function(b1,b2,b3){
  mu_f(b1,b2)*(exp(sa_f(b2))-1)/((mu_f(b1,b2)*(exp(sa_f(b2))-1))+k_f(b3)+1)
}

#'  ICC standard error for Negative Binomial with a additive extradispersion. Repeatability setting
#'
#'  ICC standard error for Negative Binomial with additive extradispersion. Repeatability setting
#' @keywords internal
#' @param b1 Intercept of GLMM
#' @param b2 Parameter related to between-subjects variance in the GLMM model
#' @param b3 Parameter related to extra-dispersion in the GLMM model
#' @param CovB Variance covariance matrix of b1, b2 and b3
#'
#' @return Scalar
#'
se_r_NB1<-function(b1,b2,b3,CovB){

  d1<-Deriv(r_NB1,"b1")
  d2<-Deriv(r_NB1,"b2")
  d3<-Deriv(r_NB1,"b3")

  dev.theta<-matrix(
    c(d1(b1,b2,b3) ,
      d2(b1,b2,b3),
      d3(b1,b2,b3))
    ,nrow=3)


  S<-array(0,c(3,3))
  S[1:3,1:3]<-CovB

  sqrt(t(dev.theta)%*%S%*%dev.theta)

}

#'  ICC for Negative Binomial with additive extradispersion. Concordance setting
#'
#'  ICC for Negative Binomial with additive extradispersion. Concordance setting
#' @keywords internal
#' @param b1 Intercept of GLMM
#' @param b2 Parameter related to between-subjects variance in the GLMM model
#' @param b3 Parameter related to extra-dispersion in the GLMM model
#' @param sb Between-methods variability
#' @return Scalar
r_NB1_2<-function(b1,b2,b3,sb){
  mu_f2(b1,b2,sb)*(exp(sa_f(b2))-1)/((mu_f2(b1,b2,sb)*(exp(sa_f(b2)+sb)-1))+k_f(b3)+1)
}

#'  ICC standard error for Negative Binomial with additive extradispersion. Repeatability setting
#'
#'  ICC standard error for Negative Binomial with additive extradispersion. Repeatability setting
#' @keywords internal
#' @param b1 Intercept of GLMM
#' @param b2 Parameter related to between-subjects variance in the GLMM model
#' @param b3 Parameter related to extra-dispersion in the GLMM model
#' @param CovB Variance covariance matrix of b1, b2 and b3
#' @param sb Between-methods variability
#' @param var.sb Standard error of sb (in variance scale).
#' @return Scalar
se_r_NB1_2<-function(b1,b2,b3,sb,CovB,var.sb){

  d1<-Deriv(r_NB1_2,"b1")
  d2<-Deriv(r_NB1_2,"b2")
  d3<-Deriv(r_NB1_2,"b3")
  d4<-Deriv(r_NB1_2,"sb")

  dev.theta<-matrix(
    c(d1(b1,b2,b3,sb) ,
      d2(b1,b2,b3,sb),
      d3(b1,b2,b3,sb),
      d4(b1,b2,b3,sb))
    ,nrow=4)


  S<-array(0,c(4,4))
  S[1:3,1:3]<-CovB
  S[4,4]<-var.sb

  sqrt(t(dev.theta)%*%S%*%dev.theta)

}

#'  ICC Zero-Inflated Poisson model. Repeatability setting
#'
#'  ICC for Zero-Inflated Poisson model. Repeatability setting
#' @keywords internal
#' @param b1 Intercept of GLMM
#' @param b2 Parameter related to between-subjects variance in the GLMM model
#' @param b3 Parameter related to extra zero probability in the GLMM model
#' @return Scalar

r_ZIP<-function(b1,b2,b3){
  pi=p_f(b3)
  (mu_f(b1,b2)*(exp(sa_f(b2))-1)*(1-pi))/(mu_f(b1,b2)*(exp(sa_f(b2))-1)+1+mu_f(b1,b2)*pi)
}

#'  ICC standard error for Zero-Inflated Poisson model. Repeatability setting
#'
#'  ICC standard error for Zero-Inflated Poisson model. Repeatability setting
#' @keywords internal
#' @param b1 Intercept of GLMM
#' @param b2 Parameter related to between-subjects variance in the GLMM model
#' @param b3 Parameter related to extra zero probability in the GLMM model
#' @param CovB Variance covariance matrix of b1,b2 and b3
#'
#' @return Scalar

se_r_ZIP<-function(b1,b2,b3,CovB){

  d1<-Deriv(r_ZIP,"b1")
  d2<-Deriv(r_ZIP,"b2")
  d3<-Deriv(r_ZIP,"b3")

  dev.theta<-matrix(
    c(d1(b1,b2,b3) ,
      d2(b1,b2,b3),
      d3(b1,b2,b3))
    ,nrow=3)

  sqrt(t(dev.theta)%*%CovB%*%dev.theta)

}

#'  ICC Zero-Inflated Poisson model. Concordance setting
#'
#'  ICC for Zero-Inflated Poisson model. Concordance setting
#' @keywords internal
#' @param b1 Intercept of GLMM
#' @param b2 Parameter related to between-subjects variance in the GLMM model
#' @param b3 Parameter related to extra zero probability in the GLMM model
#' @param sb Between-methods variability
#' @return Scalar
r_ZIP2<-function(b1,b2,b3,sb){
  pi=p_f(b3)
  (mu_f2(b1,b2,sb)*(exp(sa_f(b2))-1)*(1-pi))/(mu_f2(b1,b2,sb)*(exp(sa_f(b2)+sb)-1)+1+mu_f2(b1,b2,sb)*pi)
}


#'  ICC standard error for Zero-Inflated Poisson model. Concordance setting
#'
#'  ICC standard error for Zero-Inflated Poisson model. Concordance setting
#' @keywords internal
#' @param b1 Intercept of GLMM
#' @param b2 Parameter related to between-subjects variance in the GLMM model
#' @param b3 Parameter related to extra zero probability in the GLMM model
#' @param sb Between-methods variability
#' @param CovB Variance covariance matrix of b1, b2 and b3.
#' @param var.sb Standard error of sb (in variance scale).
#' @return Scalar

se_r_ZIP2<-function(b1,b2,b3,sb,CovB,var.sb){

  d1<-Deriv(r_ZIP2,"b1")
  d2<-Deriv(r_ZIP2,"b2")
  d3<-Deriv(r_ZIP2,"b3")
  d4<-Deriv(r_ZIP2,"sb")

  dev.theta<-matrix(
    c(d1(b1,b2,b3,sb) ,
      d2(b1,b2,b3,sb),
      d3(b1,b2,b3,sb),
      d4(b1,b2,b3,sb))
    ,nrow=4)


  S<-array(0,c(4,4))
  S[1:3,1:3]<-CovB
  S[4,4]<-var.sb

  sqrt(t(dev.theta)%*%S%*%dev.theta)

}

#'  ICC for Zero-Inflated Negative Binomial with additive extradispersion. Repeatability setting
#'
#'  ICC for Zero-Inflated Negative Binomial with additive extradispersion. Repeatability setting
#' @keywords internal
#' @param b1 Intercept of GLMM
#' @param b2 Parameter related to between-subjects variance in the GLMM model
#' @param b3 Parameter related to extra-dispersion in the GLMM model
#' @param b4 Parameter related to extra zero probability in the GLMM model
#' @return Scalar
r_ZINB1<-function(b1,b2,b3,b4){
  (mu_f(b1,b2)*(exp(sa_f(b2))-1)*(1-p_f(b4)))/(mu_f(b1,b2)*(exp(sa_f(b2))-1)+1+k_f(b3)+mu_f(b1,b2)*p_f(b4))
}

#'  ICC standard error for Zero-Inflated Negative Binomial with a additive extradispersion. Repeatability setting
#'
#'  ICC standard error for Zero-Inflated Negative Binomial with additive extradispersion. Repeatability setting
#' @keywords internal
#' @param b1 Intercept of GLMM
#' @param b2 Parameter related to between-subjects variance in the GLMM model
#' @param b3 Parameter related to extra-dispersion in the GLMM model
#' @param b4 Parameter related to extra zero probability in the GLMM model
#' @param CovB Variance covariance matrix of b1, b2, b3 and b4
#' @return Scalar
#'
se_r_ZINB1<-function(b1,b2,b3,b4,CovB){

  d1<-Deriv(r_ZINB1,"b1")
  d2<-Deriv(r_ZINB1,"b2")
  d3<-Deriv(r_ZINB1,"b3")
  d4<-Deriv(r_ZINB1,"b4")

  dev.theta<-matrix(
    c(d1(b1,b2,b3,b4) ,
      d2(b1,b2,b3,b4),
      d3(b1,b2,b3,b4),
      d4(b1,b2,b3,b4))
    ,nrow=4)


  S<-array(0,c(4,4))
  S[1:4,1:4]<-CovB

  sqrt(t(dev.theta)%*%S%*%dev.theta)

}


#'  ICC for Zero-Inflated Negative Binomial with additive extradispersion. Concordance setting
#'
#'  ICC for Zero-Inflated Negative Binomial with additive extradispersion. Concordance setting
#' @keywords internal
#' @param b1 Intercept of GLMM
#' @param b2 Parameter related to between-subjects variance in the GLMM model
#' @param b3 Parameter related to extra-dispersion in the GLMM model
#' @param b4 Parameter related to extra zero probability in the GLMM model
#' @param sb Between-methods variability
#' @return Scalar
r_ZINB1_2<-function(b1,b2,b3,b4,sb){
  (mu_f2(b1,b2,sb)*(exp(sa_f(b2))-1)*(1-p_f(b4)))/(mu_f2(b1,b2,sb)*(exp(sa_f(b2)+sb)-1)+1+k_f(b3)+mu_f2(b1,b2,sb)*p_f(b4))
}


#'  ICC standard error for Zero-Inflated Negative Binomial with additive extradispersion. Repeatability setting
#'
#'  ICC standard error for Zero-Inflated Negative Binomial with additive extradispersion. Repeatability setting
#' @keywords internal
#' @param b1 Intercept of GLMM
#' @param b2 Parameter related to between-subjects variance in the GLMM model
#' @param b3 Parameter related to extra-dispersion in the GLMM model
#' @param b4 Parameter related to extra zero probability in the GLMM model
#' @param sb Between-methods variability
#' @param CovB Variance covariance matrix of b1, b2, b3 and b4
#' @param var.sb Standard error of sb (in variance scale).
#' @return Scalar
se_r_ZINB1_2<-function(b1,b2,b3,b4,sb,CovB,var.sb){

  d1<-Deriv(r_ZINB1_2,"b1")
  d2<-Deriv(r_ZINB1_2,"b2")
  d3<-Deriv(r_ZINB1_2,"b3")
  d4<-Deriv(r_ZINB1_2,"b4")
  d5<-Deriv(r_ZINB1_2,"sb")


  dev.theta<-matrix(
    c(d1(b1,b2,b3,b4,sb) ,
      d2(b1,b2,b3,b4,sb),
      d3(b1,b2,b3,b4,sb),
      d4(b1,b2,b3,b4,sb),
      d5(b1,b2,b3,b4,sb))
    ,nrow=5)


  S<-array(0,c(5,5))
  S[1:4,1:4]<-CovB
  S[5,5]<-var.sb

  sqrt(t(dev.theta)%*%S%*%dev.theta)

}

#'  ICC for Zero-Inflated Negative Binomial with proportional extradispersion. Repeatability setting
#'
#'  ICC for Zero-Inflated Negative Binomial with proportional extradispersion. Repeatability setting
#' @keywords internal
#' @param b1 Intercept of GLMM
#' @param b2 Parameter related to between-subjects variance in the GLMM model
#' @param b3 Parameter related to extra-dispersion in the GLMM model
#' @param b4 Parameter related to extra zero probability in the GLMM model
#' @return Scalar
r_ZINB2<-function(b1,b2,b3,b4){
  K<-(1+k_f(b3))/(k_f(b3))
  (mu_f(b1,b2)*(exp(sa_f(b2))-1)*(1-p_f(b4)))/(mu_f(b1,b2)*(K*exp(sa_f(b2))-1)+1+mu_f(b1,b2)*p_f(b4))
}


#'  ICC standard error for Zero-Inflated Negative Binomial with proportional extradispersion. Repeatability setting
#'
#'  ICC standard error for Zero-Inflated Negative Binomial with proportional extradispersion. Repeatability setting
#' @keywords internal
#' @param b1 Intercept of GLMM
#' @param b2 Parameter related to between-subjects variance in the GLMM model
#' @param b3 Parameter related to extra-dispersion in the GLMM model
#' @param b4 Parameter related to extra zero probability in the GLMM model
#' @param CovB Variance covariance matrix of b1, b2, b3 and b4
#' @return Scalar

se_r_ZINB2<-function(b1,b2,b3,b4,CovB){

  d1<-Deriv(r_ZINB2,"b1")
  d2<-Deriv(r_ZINB2,"b2")
  d3<-Deriv(r_ZINB2,"b3")
  d4<-Deriv(r_ZINB2,"b4")


  dev.theta<-matrix(
    c(d1(b1,b2,b3,b4) ,
      d2(b1,b2,b3,b4),
      d3(b1,b2,b3,b4),
      d4(b1,b2,b3,b4))
    ,nrow=4)


  S<-array(0,c(4,4))
  S[1:4,1:4]<-CovB

  sqrt(t(dev.theta)%*%S%*%dev.theta)

}

#'  ICC for Zero-Inflated Negative Binomial with proportional extradispersion. Concordance setting
#'
#'  ICC for Zero-Inflated Negative Binomial with proportional extradispersion. Concordance setting
#' @keywords internal
#' @param b1 Intercept of GLMM
#' @param b2 Parameter related to between-subjects variance in the GLMM model
#' @param b3 Parameter related to extra-dispersion in the GLMM model
#' @param b4 Parameter related to extra zero probability in the GLMM model
#' @param sb Between-methods variability
#' @return Scalar

r_ZINB2_2<-function(b1,b2,b3,b4,sb){
  K<-(1+k_f(b3))/(k_f(b3))
  (mu_f2(b1,b2,sb)*(exp(sa_f(b2))-1)*(1-p_f(b4)))/(mu_f2(b1,b2,sb)*(K*exp(sa_f(b2)+sb)-1)+1+mu_f(b1,b2)*p_f(b4))

}

#'  ICC standard error for Zero-Inflated Negative Binomial with proportional extradispersion. Repeatability setting
#'
#'  ICC standard error for Zero-Inflated Negative Binomial with proportional extradispersion. Repeatability setting
#' @keywords internal
#' @param b1 Intercept of GLMM
#' @param b2 Parameter related to between-subjects variance in the GLMM model
#' @param b3 Parameter related to extra-dispersion in the GLMM model
#' @param b4 Parameter related to extra zero probability in the GLMM model
#' @param sb Between-methods variability
#' @param CovB Variance covariance matrix of b1, b2, b3 and b4
#' @param var.sb Standard error of sb (in variance scale).
#' @return Scalar

se_r_ZINB2_2<-function(b1,b2,b3,b4,sb,CovB,var.sb){

  d1<-Deriv(r_ZINB2_2,"b1")
  d2<-Deriv(r_ZINB2_2,"b2")
  d3<-Deriv(r_ZINB2_2,"b3")
  d4<-Deriv(r_ZINB2_2,"b4")
  d5<-Deriv(r_ZINB2_2,"sb")


  dev.theta<-matrix(
    c(d1(b1,b2,b3,b4,sb) ,
      d2(b1,b2,b3,b4,sb),
      d3(b1,b2,b3,b4,sb),
      d4(b1,b2,b3,b4,sb),
      d5(b1,b2,b3,b4,sb))
    ,nrow=5)


  S<-array(0,c(5,5))
  S[1:4,1:4]<-CovB
  S[5,5]<-var.sb

  sqrt(t(dev.theta)%*%S%*%dev.theta)

}


#' Intraclass correlation coefficient (ICC) for count data
#'
#' Estimates the intraclass correlation coefficient (ICC) for count data
#' @export
#' @param data A data frame containing at least two columns: outcome and subject identifier. In case of estimating the ICC for concordance setting, a third column with the method identifier will be needed.
#' @param y Character string indicating the name of the outcome column in the dataset.
#' @param id Character string indicating the name of the subjects column in the dataset.
#' @param met Character string indicating the name of the methods column in the dataset. Only needed in the concordance setting.
#' @param type Character string. It chooses the setting in which the ICC should be estimated. Valid values are: \code{"rep"} (default) for repeatability setting; \code{"con"} for concordance setting.
#' @param fam Character string. The within-subjects pdf to use. Valid options are: \code{"poisson"} (default) for Poisson pdf; \code{"nbinom1"} for Negative Binomial pdf with variance increasing linearly with the mean; \code{"nbinom2"} for Negative Binomial pdf with variance increasing quadratically with the mean; \code{"zip"} for zero-inflated Poisson pdf; \code{"zinb1"} for zero-inflated Negative Binomial nbinom1 pdf; \code{"zinb2"} for zero-inflated Negative Binomial nbinom2 pdf;
#' @param conf Confidence level for the confidence interval estimation. Default value is set to 95%.
#' @return An object of class *iccc*.The output is a list with the following components:
#' \itemize{
#'   \item *model*. An object of class glmmTMB. The estimated generalized linear mixed model.
#'   \item *ICC*. Estimate of the ICC, its standard error and confidence interval.
#'   \item *varcomp*. Variance components and parameters related to ICC expression.
#' }
#' @details The intraclass correlation coefficient (ICC) is estimated using the variance components of a generalized linear mixed model (GLMM) (Carrasco, 2010).
#'
#' The GLMM is estimated using the *glmmTMB* package (Brooks et al. 2017). The ICC standard error is estimated by applying the delta method (Ver Hoef, 2012) using the variance-covariance matrix of parameters involved in the ICC estimate.
#'
#' The parameters involved in the estimation of the ICC depends on the within-subjects pdf family chosen: the between-subjects variance, the between-methods variability (in case of concordance analysis), and parameters implicated in the within-subjects family chosen.
#' In all cases the output includes the overall expectation identified as *mu*; the between-subjects variance named as *BSVar* (log-scale); the between-methods variability (in case of concordance analysis) named as *BMVar* (log-scale).
#'
#' In the Negative Binomial pdf with variance linearly increasing with the mean (Hardin and Hilbe, 2007),
#' \deqn{Var(y_i)=\mu_i*(1+r)}
#' and Negative Binomial pdf with variance quadratically increasing with the mean (Hardin and Hilbe, 2007)
#' \deqn{Var(y_i)=\mu_i*(1+r*\mu_i)}
#' the extra-dispersion parameter *r* is included in the output.
#'
#' For zero-inflated models, the probability of observing an extra zero is included in the output as *pi*.
#'
#'
#' @references{
#' Brooks, M. E., Kristensen, K., van Benthem, K. J., Magnusson, A., Berg, C. W., Nielsen, A., Skaug, H. J., Mächler, M. and Bolker, B. M. (2017). glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling. The R Journal, 9(2), 378–400.
#'
#' Carrasco, J. (2010). A Generalized Concordance Correlation Coefficient Based on the Variance Components Generalized Linear Mixed Models for Overdispersed Count Data. Biometrics, 66(3), 897-904.
#'
#' W. Hardin and J. Hilbe. (2007). Generalized Linear Models and Extensions. Stata Press.
#'
#' Ver Hoef, J.M. (2012) Who Invented the Delta Method?, The American Statistician, 66:2, 124-127,
#' }
#' @examples
#' # Poisson model. Repeatability setting.
#' iccpois<-icc_counts(Grimso,y="Tot",id="TransectID")
#' # Negative Binomial with proportional extra-dispersion. Concordance setting
#' iccnb2<-icc_counts(AF,y="y",id="id",met="met",type="con",fam="nbinom2")
#' # Zero-inflated Poisson model. Repeatability setting
#' icczip<-icc_counts(EPP,y="Social",id="id",fam="zip")


icc_counts<-function(data,y,id,met=NULL,type=c("rep","con"),
                     fam=c("poisson","nbinom1","nbinom2","zip","zinb1","zinb2"),
                     conf=0.95){

  type<-match.arg(type,choices=c("rep","con"))
  fam<-match.arg(fam,choices=c("poisson","nbinom1","nbinom2","zip","zinb1","zinb2"))
  x<-get_data(data,y=y,id=id,met=met)


  if (type=="con"){

    if (fam=="poisson"){
      m1<-fit_model_con(x)
      est1<-r_est_con(x,m1)
      varcomp=est1[3:5]
    }

    if (fam=="nbinom2"){
      m1<-fit_model_con(x,fam="nbinom2")
      est1<-r_est_con(x,m1)
      varcomp=est1[3:6]
    }
    if (fam=="nbinom1"){
      m1<-fit_model_con(x,fam="nbinom1")
      est1<-r_est_con(x,m1)
      varcomp=est1[3:6]
    }

    if (fam=="zip"){
      m1<-fit_model_con(x,fam="zip")
      est1<-r_est_con_zi(x,m1)
      varcomp=est1[3:6]

    }

    if (fam=="zinb1"){
      m1<-fit_model_con(x,fam="zinb1")
      est1<-r_est_con_zi(x,m1)
      varcomp=est1[3:7]

    }

    if (fam=="zinb2"){
      m1<-fit_model_con(x,fam="zinb2")
      est1<-r_est_con_zi(x,m1)
      varcomp=est1[3:7]

    }


  }

  if (type=="rep"){

    if (fam=="poisson"){
      m1<-fit_model_rep(x)
      est1<-r_est_rep(x,m1)
      varcomp=est1[3:4]

    }

    if (fam=="nbinom2"){
      m1<-fit_model_rep(x,fam="nbinom2")
      est1<-r_est_rep(x,m1)
      varcomp=est1[3:5]
    }

    if (fam=="nbinom1"){
      m1<-fit_model_rep(x,fam="nbinom1")
      est1<-r_est_rep(x,m1)
      varcomp=est1[3:5]
    }

    if (fam=="zip"){
      m1<-fit_model_rep(x,fam="zip")
      est1<-r_est_rep_zi(x,m1)
      varcomp=est1[3:5]

    }

    if (fam=="zinb1"){
      m1<-fit_model_rep(x,fam="zinb1")
      est1<-r_est_rep_zi(x,m1)
      varcomp=est1[3:6]

    }

    if (fam=="zinb2"){
      m1<-fit_model_rep(x,fam="zinb2")
      est1<-r_est_rep_zi(x,m1)
      varcomp=est1[3:6]

    }

  }


  rownames(varcomp)<-""

  out<-list(model=m1,ICC=ic_r(est1$ICC,est1$`SE ICC`^2,conf=conf),
            varcomp=varcomp)
  class(out)<-"iccc"
  return(out)
}



#' Bland-Altman plot
#'
#' Draws the Bland-Altman plot. The differences among pair of data from the same subject
#'  is represented on y-axis. The mean of data from the same subject is represented
#'  on x-axis. Additionally, a bar plot with the proportions of differences
#'  can be drawn.
#'
#' @export
#' @param data A data frame containing at least two columns: outcome and subject identifier.
#' @param y Character string indicating the name of the outcome column in the data set.
#' @param id Character string indicating the name of the subjects column in the data set.
#' @param rm Optional. Character string indicating the name of  column that stands for the repeated measurements from the same subjects in the dataset.
#' Only needed to identify the differences in the Bland-Altman plot.
#' @param type Character. Which plot has to be drawn? Default option is Bland-Altman plot ("BA" option). Alternatively, the bar plot of the proportion of the differences can be created ("bars" option).
#' @return A list with the following components:
#' \itemize{
#'   \item *plot*. An object of class ggplot. The plot generated.
#'   \item *data*. An object of class dataframe that contains the data used to generated the plot.
#' }
#' @examples
#' plot_BA(EPP,y="Social",id="id")
#' plot_BA(EPP,y="Social",id="id",rm="Year")
#' plot_BA(EPP,y="Social",id="id",type="bars")
#
#'
plot_BA<-function(data,y,id,rm=NULL,type=c("BA","bars")){

  type<-match.arg(type)
  x<-get_data_plot(data,y=y,id=id,rm=rm)

  if (is.null(rm)==T) {
    aux<-x %>% group_by(id) %>% dplyr::reframe(Diff = combn(y,diff,m=2),
                                          m = mean(y))

    my<-max(abs(aux$Diff))
    mx<-max(aux$m)
    mnx<-min(abs(aux$m))
    rx<-mx-mnx

    g1<-ggplot(data=aux,aes(x=m,y=(Diff))) + geom_point() +
      xlim(mnx-0.05*rx,mx+0.05*rx)+
      ylim(my*(-1),my) + geom_hline(yintercept = 0) +
      xlab("Mean") + ylab("Difference")

  }


  if ( (is.null(rm)==F) &  (is.null(x$rm)==F)){

    aux<-x %>% group_by(id) %>% dplyr::reframe(Diff = combn(y,diff,m=2),
                                          m = mean(y),
                                          lev1=combn(rm,m=2, function(x) x[1]),
                                          lev2=combn(rm,m=2, function(x) x[2])) %>%
      mutate(Levels=paste(lev1,"-",lev2))

    my<-max(abs(aux$Diff))
    mx<-max(aux$m)
    mnx<-min(abs(aux$m))
    rx<-mx-mnx

    g1<-ggplot(data=aux,aes(x=m,y=Diff,color=Levels)) + geom_point() +
      xlim(mnx-0.05*rx,mx+0.05*rx)+
      ylim(my*(-1),my) + geom_hline(yintercept = 0) +
      xlab("Mean") + ylab("Difference")

  }

  if ( (is.null(rm)==F) &  (is.null(x$rm)==T)) stop("Incorrect name for repeated measures variable")

  if ( (is.null(rm)==T) & (type=="bars") ) {

    aux2<-aux %>% group_by(Diff) %>% dplyr::reframe(n=n(),per=n/nrow(aux))
    g1<-ggplot(data=aux2,aes(x=Diff,y=per)) + geom_bar(stat="identity") +
      xlab("Difference") + ylab("Proportion")

  }


    if ( (is.null(rm)==F) & (type=="bars") ) stop("Repeated measures variable option only available for BA plot")

  out<-list(plot=g1,data=aux)
  print(g1)
  return(invisible(out))
}




