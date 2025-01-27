% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auxiliar_functions.R
\name{VarComp}
\alias{VarComp}
\title{GLMM variance components}
\usage{
VarComp(x, digits = getOption("digits"))
}
\arguments{
\item{x}{An object of class \emph{iccc}}

\item{digits}{Number of digits to print}
}
\value{
A vector with the variance components and related parameters
}
\description{
Prints the GLMM variance components and related parameters to estimate the intraclass correlation coefficient (ICC)
}
\examples{
# Poisson model. Repeatability setting.
iccpois<-icc_counts(Grimso,y="Tot",id="TransectID")
VarComp(iccpois)
# Negative Binomial with proportional extra-dispersion. Concordance setting
iccnb2<-icc_counts(AF,y="y",id="id",met="met",type="con",fam="nbinom2")
VarComp(iccnb2)
# Zero-inflated Poisson model. Repeatability setting
icczip<-icc_counts(EPP,y="Social",id="id",fam="zip")
VarComp(icczip)
}
\seealso{
\code{\link[=icc_counts]{icc_counts()}}
}
