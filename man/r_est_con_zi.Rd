% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_functions.R
\name{r_est_con_zi}
\alias{r_est_con_zi}
\title{ICC concordance for zero-inflated models}
\usage{
r_est_con_zi(data, model)
}
\arguments{
\item{data}{A data frame containing at least three columns: outcome (named as y), subject identifier (named as id) and method identifier (named as met).}

\item{model}{A glmmTMB obejct.}
}
\value{
A data frame with ICC estimate, its standard error and variance components.
}
\description{
Estimates the ICC for concordance setting for zero-inflated models (poisson, nbinom1 and nbinom2 families).
}
\keyword{internal}
