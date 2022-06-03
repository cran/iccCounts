#' @importFrom stats poisson pnbinom ppois predict qnorm runif simulate quantile
#' @importFrom VGAM pzipois pzinegbin
#' @importFrom gridExtra grid.arrange
#' @importFrom Deriv Deriv
#' @importFrom utils combn globalVariables
#' @import glmmTMB dplyr ggplot2


utils::globalVariables(c("m", "Diff", "lev1","lev2","Levels","per"))

# Transformations

f1<-function(s) exp(s)
f2<-function(s) exp(s)/(1+exp(s))
f3<-function(s) exp(s)^2
f4<-function(b1,b2) log(mu_f(b1,b2))

# Expectations and other parameters
mu_f<-function(b1,b2) f1(b1+f3(b2)/2)
mu_f2<-function(b1,b2,sb) f1(b1+(f3(b2)+sb)/2)

sa_f<-function(b2) f3(b2)
k_f<-function(b3) f1(b3) # gof parameter negbinom
p_f<-function(b4) f2(b4) # inflation parameter ZI



# Builds L matrix for 2x2 comparison of means
Lmat<-function(k){
  aux1<-diag(k)*0
  aux1[1,1:k]<-1

  for (i in 1:(k-1)) aux1[1+i,1+i]=1

  L=array(0,dim=c(k,k*(k-1)/2))
  cont=0
  for (i in 1:k){
    for (j in 1:(k-1)){
      if (i>j){
        cont=cont+1
        L[,cont]=aux1[,i]-aux1[,j]
      }
    }
  }
  return(L)
}



get_data<-function(data,y,id,met=NULL){

  data<-as.data.frame(data)
  vid<-data[,id]
  vy<-data[,y]

  if (is.null(met)==F){
    vmet<-data[,met]
    vmet<-as.factor(vmet)
    x<-data.frame(vy,vid,vmet)
    names(x)<-c("y","id","met")
    if (length(unique(x$met))<2) stop("Less than two methods")
    return(x)
  }  else {
    x<-data.frame(vy,vid)
    names(x)<-c("y","id")
    return(x)
  }

}


get_data_plot<-function(data,y,id,rm=NULL){

  data<-as.data.frame(data)
  vid<-data[,id]
  vy<-data[,y]

  if (is.null(rm)==F){
    vrm<-data[,rm]
    vrm<-as.factor(vrm)
    x<-data.frame(vy,vid,vrm)
    names(x)<-c("y","id","rm")
    if (length(unique(x$rm))<2) stop("Less than two repeated measures")
    return(x)
  }  else {
    x<-data.frame(vy,vid)
    names(x)<-c("y","id")
    return(x)
  }

}


sb_f<-function(k,b,varB){

  L<-Lmat(k)   # Building L matrix

  difmed=t(L)%*%b  # Calculating observers sum of squares

  A<-L%*%t(L)
  aux2<-(t(difmed)%*%difmed)-sum(diag(A%*%varB))
  sb<-max(aux2/(k*(k-1)),0)
  return(sb)
}


# Variance of sb
var_sb<-function(k,b,varB){
  L<-Lmat(k)    # Building L matrix
  A<-L%*%t(L)
  xx<-( 2*sum(diag( ( (A%*%varB)**2 ) ) )+ 4*( t(b)%*%A%*%varB%*%A%*%b ) ) /((k*(k-1))**2)
  return(xx)
}


ic_r<-function(r,var.r,conf=0.95){

  alpha<-1-conf
  z<-0.5*log((1+r)/(1-r))
  var.z<-var.r/(((1+r)**2)*((1-r)**2))
  ic.z<-z+c(-1,1)*as.vector(qnorm(1-(alpha/2))*sqrt(var.z))
  ic.r<-(exp(2*ic.z)-1)/(exp(2*ic.z)+1)
  out<-matrix(c(r,sqrt(var.r),ic.r),nrow=1)

  colnames(out)<-c("ICC","SE ICC",paste(conf*100,"% CI LL",sep=""),
                   paste(conf*100,"% CI UL",sep=""))
  return(out)

}


varB<-function(model,k) model$sdr$cov.fixed[1:k,1:k]


#' Prints the intraclass correlation coefficient
#'
#' Prints the intraclass correlation coefficient (ICC)
#'
#' @param x An object of class *iccc*
#' @param digits Number of digits to print
#' @return A vector with the ICC estimate, its standard error and confidence interval.
#' @export
#' @examples
#' # Poisson model. Repeatability setting.
#' iccpois<-icc_counts(Grimso,y="Tot",id="TransectID")
#' ICC(iccpois)
#' # Negative Binomial with proportional extra-dispersion. Concordance setting
#' iccnb2<-icc_counts(AF,y="y",id="id",met="met",type="con",fam="nbinom2")
#' ICC(iccnb2)
#' # Zero-inflated Poisson model. Repeatability setting
#' icczip<-icc_counts(EPP,y="Social",id="id",fam="zip")
#' ICC(icczip)
#' @seealso [icc_counts()]
ICC<-function(x,digits = getOption("digits")){

  fam<-x$model$modelInfo$family$family
  zi<-x$model$modelInfo$allForm$ziformula
  if (zi==~1) ot<-"zero inflated" else ot<-""

  message("Model: ",paste(fam,ot))

  return(round(x$ICC,digits = digits))

}

#' GLMM variance components
#'
#' Prints the GLMM variance components and related parameters to estimate the intraclass correlation coefficient (ICC)
#'
#' @param x An object of class *iccc*
#' @param digits Number of digits to print
#' @return A vector with the variance components and related parameters
#' @export
#' @examples
#' # Poisson model. Repeatability setting.
#' iccpois<-icc_counts(Grimso,y="Tot",id="TransectID")
#' VarComp(iccpois)
#' # Negative Binomial with proportional extra-dispersion. Concordance setting
#' iccnb2<-icc_counts(AF,y="y",id="id",met="met",type="con",fam="nbinom2")
#' VarComp(iccnb2)
#' # Zero-inflated Poisson model. Repeatability setting
#' icczip<-icc_counts(EPP,y="Social",id="id",fam="zip")
#' VarComp(icczip)
#' @seealso [icc_counts()]
VarComp<-function(x,digits = getOption("digits")){

  fam<-x$model$modelInfo$family$family
  zi<-x$model$modelInfo$allForm$ziformula
  if (zi==~1) ot<-"zero inflated" else ot<-""
  message("Model: ",paste(fam,ot))
  return(round(x$varcomp,digits = digits))

}



count_zero<-function(x) sum(x==0)

#' Goodness of fit plots
#'
#' Draws the plots to assess the goodness of fit
#'
#' @param x An object of class *GOF* generated by *GOF_check* function.
#' @param type Which plot to draw. Values: *all* (default); *envelope* for envelopes of randomized quantile residuals; *dispersion* for plot to assess the dispersion; *zeros* for plot to assess the zero inflation.
#' @param ... Ignore
#' @method plot GOF
#' @export
#' @examples
#' \donttest{
#' # Poisson model. Repeatability setting.
#' iccpois<-icc_counts(EPP,y="Social",id="id")
#' iccpois.gof<-GOF_check(iccpois)
#' plot(iccpois.gof)
#' plot(iccpois.gof,type="envelope")
#' plot(iccpois.gof,type="dispersion")
#' plot(iccpois.gof,type="zeros")
#' }
#' @seealso [GOF_check()]
plot.GOF<-function(x,type=c("all","envelope","dispersion","zeros"),...){
  type<-match.arg(type,choices=c("all","envelope","dispersion","zeros"))

  if (type=="all") grid.arrange(x$plot_env,x$plot_var,x$plot_zi)
  if (type=="envelope") return(x$plot_env)
  if (type=="dispersion") return(x$plot_var)
  if (type=="zeros") return(x$plot_zi)

}



