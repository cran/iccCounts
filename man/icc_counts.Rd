% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_functions.R
\name{icc_counts}
\alias{icc_counts}
\title{Intraclass correlation coefficient (ICC) for count data}
\usage{
icc_counts(
  data,
  y,
  id,
  met = NULL,
  type = c("rep", "con"),
  fam = c("poisson", "nbinom1", "nbinom2", "zip", "zinb1", "zinb2"),
  conf = 0.95
)
}
\arguments{
\item{data}{A data frame containing at least two columns: outcome and subject identifier. In case of estimating the ICC for concordance setting, a third column with the method identifier will be needed.}

\item{y}{Character string indicating the name of the outcome column in the dataset.}

\item{id}{Character string indicating the name of the subjects column in the dataset.}

\item{met}{Character string indicating the name of the methods column in the dataset. Only needed in the concordance setting.}

\item{type}{Character string. It chooses the setting in which the ICC should be estimated. Valid values are: \code{"rep"} (default) for repeatability setting; \code{"con"} for concordance setting.}

\item{fam}{Character string. The within-subjects pdf to use. Valid options are: \code{"poisson"} (default) for Poisson pdf; \code{"nbinom1"} for Negative Binomial pdf with variance increasing linearly with the mean; \code{"nbinom2"} for Negative Binomial pdf with variance increasing quadratically with the mean; \code{"zip"} for zero-inflated Poisson pdf; \code{"zinb1"} for zero-inflated Negative Binomial nbinom1 pdf; \code{"zinb2"} for zero-inflated Negative Binomial nbinom2 pdf;}

\item{conf}{Confidence level for the confidence interval estimation. Default value is set to 95\%.}
}
\value{
An object of class \emph{iccc}.The output is a list with the following components:
\itemize{
\item \emph{model}. An object of class glmmTMB. The estimated generalized linear mixed model.
\item \emph{ICC}. Estimate of the ICC, its standard error and confidence interval.
\item \emph{varcomp}. Variance components and parameters related to ICC expression.
}
}
\description{
Estimates the intraclass correlation coefficient (ICC) for count data
}
\details{
The intraclass correlation coefficient (ICC) is estimated using the variance components of a generalized linear mixed model (GLMM) (Carrasco, 2010).

The GLMM is estimated using the \emph{glmmTMB} package (Brooks et al. 2017). The ICC standard error is estimated by applying the delta method (Ver Hoef, 2012) using the variance-covariance matrix of parameters involved in the ICC estimate.

The parameters involved in the estimation of the ICC depends on the within-subjects pdf family chosen: the between-subjects variance, the between-methods variability (in case of concordance analysis), and parameters implicated in the within-subjects family chosen.
In all cases the output includes the overall expectation identified as \emph{mu}; the between-subjects variance named as \emph{BSVar} (log-scale); the between-methods variability (in case of concordance analysis) named as \emph{BMVar} (log-scale).

In the Negative Binomial pdf with variance linearly increasing with the mean (Hardin and Hilbe, 2007),
\deqn{Var(y_i)=\mu_i*(1+r)}
and Negative Binomial pdf with variance quadratically increasing with the mean (Hardin and Hilbe, 2007)
\deqn{Var(y_i)=\mu_i*(1+r*\mu_i)}
the extra-dispersion parameter \emph{r} is included in the output.

For zero-inflated models, the probability of observing an extra zero is included in the output as \emph{pi}.
}
\examples{
# Poisson model. Repeatability setting.
iccpois<-icc_counts(Grimso,y="Tot",id="TransectID")
# Negative Binomial with proportional extra-dispersion. Concordance setting
iccnb2<-icc_counts(AF,y="y",id="id",met="met",type="con",fam="nbinom2")
# Zero-inflated Poisson model. Repeatability setting
icczip<-icc_counts(EPP,y="Social",id="id",fam="zip")
}
\references{
{
Brooks, M. E., Kristensen, K., van Benthem, K. J., Magnusson, A., Berg, C. W., Nielsen, A., Skaug, H. J., Mächler, M. and Bolker, B. M. (2017). glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling. The R Journal, 9(2), 378–400.

Carrasco, J. (2010). A Generalized Concordance Correlation Coefficient Based on the Variance Components Generalized Linear Mixed Models for Overdispersed Count Data. Biometrics, 66(3), 897-904.

W. Hardin and J. Hilbe. (2007). Generalized Linear Models and Extensions. Stata Press.

Ver Hoef, J.M. (2012) Who Invented the Delta Method?, The American Statistician, 66:2, 124-127,
}
}
