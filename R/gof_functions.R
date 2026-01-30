

Pres_ZIP<-function(y,mu,p)  (y-mu)/sqrt(mu*(1+p*mu/(1-p)))
Pres_ZINB1<-function(y,mu,p,r)  (y-mu)/sqrt(mu*(1+r+p*mu/(1-p)))
Pres_ZINB2<-function(y,mu,p,r) (y-mu)/sqrt(mu+(r*(mu^2)/(1-p))+((mu/(1-p))^2)*(p^2+p))


#' RQR for Poisson GLMM
#' Computes the randomized quantile residuals for Poisson GLMM
#'
#' @keywords internal
#'
#' @param x An object of clas *iccc*.
#' @return A vector with the residuals.
rqr_pois<-function(x){
  mu<-predict(x$model,type="response")
  y<-x$model$frame$y
  rqr<-qnorm(runif(length(y),ppois(y-1,mu),ppois(y,mu)))
  rqr<-ifelse(rqr=="Inf",5,rqr)
  rqr<-ifelse(rqr=="-Inf",-5,rqr)
  return(rqr)
}


#' RQR for NegBin1 GLMM
#' Computes the randomized quantile residuals for NegBin1 GLMM
#'
#' @keywords internal
#' @param x An object of clas *iccc*.
#' @return A vector with the residuals.
rqr_nb1<-function(x){
  mu<-predict(x$model,type="response")
  y<-x$model$frame$y
  x$varcomp
  r<-as.numeric(x$varcomp[3])
  rqr<-qnorm(runif(length(y),pnbinom(y-1,mu=mu,size=mu/r),pnbinom(y,mu=mu,size=mu/r)))
  rqr<-ifelse(rqr=="Inf",5,rqr)
  rqr<-ifelse(rqr=="-Inf",-5,rqr)
  return(rqr)
}

#' RQR for NegBin2 GLMM
#' Computes the randomized quantile residuals for NegBin2 GLMM
#'
#' @keywords internal
#' @param x An object of clas *iccc*.
#' @return A vector with the residuals.
rqr_nb2<-function(x){
  mu<-predict(x$model,type="response")
  y<-x$model$frame$y
  r<-as.numeric(x$varcomp[3])
  rqr<-qnorm(runif(length(y),pnbinom(y-1,mu=mu,size=1/r),pnbinom(y,mu=mu,size=1/r)))
  rqr<-ifelse(rqr=="Inf",5,rqr)
  rqr<-ifelse(rqr=="-Inf",-5,rqr)
  return(rqr)
}

#' RQR for ZIP GLMM
#' Computes the randomized quantile residuals for ZIP GLMM
#' @keywords internal
#' @param x An object of clas *iccc*.
#' @return A vector with the residuals.
rqr_zip<-function(x){
  mu<-predict(x$model,type="response")
  y<-x$model$frame$y
  p<-as.numeric(x$varcomp[3])
  rqr<-qnorm(runif(length(y),pzipois(y-1,lambda=mu,pstr0=p),pzipois(y,lambda=mu,pstr0=p)))
  rqr<-ifelse(rqr=="Inf",5,rqr)
  rqr<-ifelse(rqr=="-Inf",-5,rqr)
  return(rqr)
}

#' RQR for ZINB1 GLMM
#' Computes the randomized quantile residuals for ZINB1 GLMM
#' @keywords internal
#' @param x An object of clas *iccc*.
#' @return A vector with the residuals.
rqr_zinb1<-function(x){
  mu<-predict(x$model,type="response")
  y<-x$model$frame$y
  r<-as.numeric(x$varcomp[3])
  p<-as.numeric(x$varcomp[4])
  rqr<-qnorm(runif(length(y),pzinegbin(y-1,munb=mu,size=mu/r),pzinegbin(y,munb=mu,size=mu/r)))
  rqr<-ifelse(rqr=="Inf",5,rqr)
  rqr<-ifelse(rqr=="-Inf",-5,rqr)
  return(rqr)
}

#' RQR for ZINB2 GLMM
#' Computes the randomized quantile residuals for ZINB2 GLMM
#' @keywords internal
#' @param x An object of clas *iccc*.
#' @return A vector with the residuals.

rqr_zinb2<-function(x){
  mu<-predict(x$model,type="response")
  y<-x$model$frame$y
  r<-as.numeric(x$varcomp[3])
  p<-as.numeric(x$varcomp[4])
  rqr<-qnorm(runif(length(y),pzinegbin(y-1,munb=mu,size=1/r),pzinegbin(y,munb=mu,size=1/r)))
  rqr<-ifelse(rqr=="Inf",5,rqr)
  rqr<-ifelse(rqr=="-Inf",-5,rqr)
  return(rqr)
}


#' RQR for GLMM
#' Computes the randomized quantile residuals for GLMM
#' @keywords internal
#' @param x An object of clas *iccc*.
#' @details Randomized quantile residuals (RQR) are computed for GLMMM. The within-cluster families considered
#' are Poisson, Negative Binomial with additive and proportional extra-dispersion and their zero-inflated extensions.
#' For further details on RQR see Dunn and Smyth (1996) and Feng et al (2020)
#' @return A vector with the residuals.
#' @references {
#'
#' Dunn PK, Smyth GK. (1996). Randomized quantile residuals. J Comput Graph Stat.5(3):236–44.
#'
#' Feng et al. (2020). A comparison of residual diagnosis tools for diagnosing regression models for count data. BMC Medical Research Methodology 20:175
#' }
rqr<-function(x){
  fam<-x$model$modelInfo$family$family

  zi<-x$model$modelInfo$allForm$ziformula

  if (zi!=~1){
    if (fam=="poisson") res<-rqr_pois(x)
    if (fam=="nbinom1") res<-rqr_nb1(x)
    if (fam=="nbinom2") res<-rqr_nb2(x)
  }

  if (zi==~1){
    if (fam=="poisson") res<-rqr_zip(x)
    if (fam=="nbinom1") res<-rqr_zinb1(x)
    if (fam=="nbinom2") res<-rqr_zinb2(x)
  }

    return(res)
}


var_res <- function(x) {
  n<-nrow(x$model$frame)
  ns<-length(unique(x$model$frame$id))
  p<-length(x$model$fit$par)
  rdf<-n-ns-p
  res<-rqr(x)
  var.res<-sum(res^2)/rdf
  out<-matrix(var.res,nrow=1)
  colnames(out)<-c("Residual Variance")
  return(out)
}

var_res_2 <- function(x,df) {
  var.res<-sum(x^2)/df
  out<-matrix(var.res,nrow=1)
  colnames(out)<-c("Residual Variance")
  return(out)
}


getEnvelope <- function(model.object,nsim=100){

  zi<-model.object$modelInfo$allForm$ziformula
  if (zi!=~1)  fam<-model.object$modelInfo$family$family
  if (zi==~1){
    fam2<-model.object$modelInfo$family$family
    if (fam2 == "poisson") fam<-"zip"
    if (fam2 == "nbinom1") fam<-"zinb1"
    if (fam2 == "nbinom2") fam<-"zinb2"
  }


  if (is.null(model.object$frame$met) ) {
    sim.data<-cbind(simulate(model.object,nsim=nsim), id=model.object$frame$id)

    message("Simulating...","\n")

    resid.matrix<-sapply(1:nsim,function(i){

      if (i %in% seq(10,nsim,10))  message(paste(c(i,"..."),sep=""),appendLF=F)
      sim.data$y<-sim.data[,i]
      est<-invisible(icc_counts(sim.data,y="y",id="id",type=c("rep"),fam=fam))
      sort(rqr(est))
    }
    )

  }

  if (is.null(model.object$frame$met)==F ) {
    sim.data<-cbind(simulate(model.object,nsim=nsim), id=model.object$frame$id, met=model.object$frame$met)

    message("Simulating...","\n")

    resid.matrix<-sapply(1:nsim,function(i){
      if (i %in% seq(10,nsim,10))  message(paste(c(i,"..."),sep=""),appendLF=F)
      sim.data$y<-sim.data[,i]
      est<-invisible(icc_counts(sim.data,y="y",id="id",met="met",type=c("con"),
                                fam=fam))
      sort(rqr(est))
    }
    )

  }


  return(list(simdata=sim.data,resid.matrix=resid.matrix))

}



getQnorm <- function(y){

  n  <- length(y)
  qu <- rep(NA, n)

  for(i in 1:n){

    qu[i] <- sum(y <= y[i])/n

  }

  qh <- qnorm(p = qu)

  return(qh)

}


#' Goodness of fit for GLMM
#'
#' Assessment of goodness of fit for GLMM
#'
#' @param x An object of clas *iccc*.
#' @param nsim Number of simulations to run. Default is set to 100.
#' @param alpha Level of significance
#' @details Randomized quantile residuals are computed for the fitted model. Simulations based on the fitted model are generated
#' and the model is refitted to each simulated dataset. Envelopes for RQR are built as the appropriate quantile (in relation to the level fo significance) of RQR from
#' the refitted models. Additionally, a test for dispersion and zero inflation are carried out by comparing the RQR dispersion and the
#' number of zeros from the original model and data to those from the refitted models and simulated data.
#' @export
#' @return An object of class *GOF* for which method *plot* is available. A list with the following components:
#' \itemize{
#' \item *plot_env*. Plot of RQR envelopes with the original RQR.
#' \item *plot_var*. Plot of the simulated RQR dispersion.
#' \item *plot_zi*. Plot of the count of zeros in the simulated datasets.
#' \item *res_var*. Dispersion of RQR from the original sample.
#' \item *pval_var*. Proportion of simulated RQR dispersion that are greater than the original dispersion that can be interpreted as a simulated P-value to check the goodness of fit on dispersion.
#' \item *zero_count*. Count of zeros in the original sample.
#' \item *pval_zi*. Proportion of simulated zero count that are greater than that of the original sample. It can be interpreted as a simulated P-value to check the hypothesis of zero-inflation.
#' }
#' @seealso [plot.GOF()], [DispersionTest()],[ZeroTest()]
#'
#' @examples
#' \donttest{
#' # Poisson model. Repeatability setting.
#' iccpois<-icc_counts(EPP,y="Social",id="id")
#' GOF_check(iccpois)
#' # Zero-inflated Poisson model. Repeatability setting
#' icczip<-icc_counts(EPP,y="Social",id="id",fam="zip")
#' GOF_check(icczip)
#' }
GOF_check <- function(x, nsim = 100, alpha = 0.05) {
  Freq <- NULL
  sim_out <- getEnvelope(x$model, nsim = nsim)

  mat <- sim_out$resid.matrix

  envelope.df <- data.frame(
    min = apply(mat, 1, FUN = quantile, probs = alpha, type = 6),
    max = apply(mat, 1, FUN = quantile, probs = 1 - alpha, type = 6),
    mean = apply(mat, 1, FUN = mean)
  )

  res <- sort(rqr(x))
  theo <- getQnorm(res)
  env <- as.data.frame(cbind(envelope.df, res, theo))

  zi <- x$model$modelInfo$allForm$ziformula
  if (zi != ~1) plot_tit <- x$model$modelInfo$family$family
  if (zi == ~1) plot_tit <- paste("Zero Inflated", x$model$modelInfo$family$family)

  g1 <- ggplot(data = env, aes(x = theo, y = res)) +
    geom_point(shape = "+") +
    geom_line(aes(x = theo, y = mean), color = "red", lty = 2) +
    geom_line(aes(x = theo, y = min), color = "black") +
    geom_line(aes(x = theo, y = max), color = "black") +
    theme_bw() +
    labs(x = "Theoretical quantiles",
         y = "Sample quantiles",
         title = paste(plot_tit, "model"))

  # Dispersion test

  n <- nrow(x$model$frame)
  ns <- length(unique(x$model$frame$id))
  p <- length(x$model$fit$par)
  rdf <- n - ns - p
  act_var <- as.numeric(var_res(x))
  out_var_res <- data.frame(theo = apply(sim_out$resid.matrix, 2, var_res_2, df = rdf))
  pval_var <- (sum(act_var < out_var_res$theo) + 1) / (nrow(out_var_res) + 1)

  # Preparem l'etiqueta fora del ggplot
  lab.g2 <- paste("Sample dispersion \n", round(act_var, 2))

  # MODIFICAT: Eliminat 'label' de l'aes() principal
  g2 <- ggplot(out_var_res, aes(x = theo)) +
    geom_density() +
    labs(y = "Density", x = "Simulated Residual Dispersion",
         title = paste(plot_tit, "model")) +
    theme_bw()

  limy <- layer_scales(g2)$y$range$range
  limx <- layer_scales(g2)$x$range$range

  # MODIFICAT: Ús de annotate() en lloc de geom_label()
  g2 <- g2 + annotate("label",
                      x = limx[1] + 0.2 * diff(limx),
                      y = limy[1] + 0.05 * diff(limy),
                      label = lab.g2,
                      size = 3)

  # ZI test
  obs_c <- count_zero(x$model$frame$y)
  y_sim <- sim_out$simdata %>% select(starts_with("sim"))
  c0 <- apply(y_sim, 2, count_zero)
  pval_zi <- (sum(obs_c <= c0) + 1) / (length(c0) + 1)
  c0 <- factor(c0, levels = min(c0):max(c0))
  z <- prop.table(table(c0))

  lab.g3 <- paste("Zeros in sample:", obs_c)

  # MODIFICAT: Eliminat 'label' de l'aes() principal
  g3 <- ggplot(as.data.frame(z), aes(c0, Freq)) +
    geom_bar(stat = "identity") +
    xlab("Simulated Zero Count") +
    labs(y = "Proportion",
         title = paste(plot_tit, "model")) +
    theme_bw()

  # MODIFICAT: Ús de annotate() en lloc de geom_label()
  g3 <- g3 + annotate("label",
                      x = levels(c0)[1],
                      y = 0.1,
                      label = lab.g3,
                      size = 3,
                      hjust = -0.5)

  out <- list(plot_env = g1, plot_var = g2, res_var = act_var, pval_var = pval_var,
              plot_zi = g3, zero_count = obs_c, pval_zi = pval_zi)

  class(out) <- "GOF"
  return(out)
}


#' Dispersion test for GLMM
#'
#' @param x An object of class *GOF* generated by *GOF_check* function.
#' @details The function prints the dispersion of sample randomized quantile residuals (RQR) and the simulated P-value.
#' @export
#' @return A vector with the sample RQR dispersion and the P-value.
#'
#' @examples
#' \donttest{
#' # Poisson model. Repeatability setting.
#' iccpois<-icc_counts(EPP,y="Social",id="id")
#' iccpois.gof<-GOF_check(iccpois)
#' DispersionTest(iccpois.gof)
#' }
#' @seealso [GOF_check()]
DispersionTest<-function(x){
  out<-data.frame(S=x$res_var,P_value=x$pval_var,row.names="")
  return(out)
}

#' Zero-Inflation test for GLMM
#'
#' @param x An object of class *GOF* generated by *GOF_check* function.
#' @details The function prints the count of zeros in the sample and the simulated P-value.
#' @return A vector with the zero count and the P-value.
#' @export
#' @examples
#' \donttest{
#' # Poisson model. Repeatability setting.
#' iccpois<-icc_counts(EPP,y="Social",id="id")
#' iccpois.gof<-GOF_check(iccpois)
#' ZeroTest(iccpois.gof)
#' # Zero-inflated Poisson model. Repeatability setting
#' icczip<-icc_counts(EPP,y="Social",id="id",fam="zip")
#' icczip.gof<-GOF_check(icczip)
#' ZeroTest(icczip.gof)
#' }
#' @seealso [GOF_check()]
ZeroTest<-function(x){
  out<-data.frame(Count=x$zero_count,P_value=x$pval_zi,row.names="")
  return(out)
}


#' @export
print.GOF <- function(x, ...) {
  cat("\nResults of Goodness of Fit (GOF) Check\n")
  cat("--------------------------------------\n")

  # Mostrem la variància
  cat(sprintf("Residual Variance: %.4f (P-value: %.4f)\n",
              x$res_var, x$pval_var))

  # Mostrem els zeros
  cat(sprintf("Zero Count:        %d     (P-value: %.4f)\n",
              x$zero_count, x$pval_zi))

  cat("--------------------------------------\n")
  cat("Plots available in the object: $plot_env, $plot_var, $plot_zi\n")

  # Retornem l'objecte de manera invisible per no tornar-lo a imprimir
  invisible(x)
}
