% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_functions.R
\name{r_ZINB2_2}
\alias{r_ZINB2_2}
\title{ICC for Zero-Inflated Negative Binomial with proportional extradispersion. Concordance setting}
\usage{
r_ZINB2_2(b1, b2, b3, b4, sb)
}
\arguments{
\item{b1}{Intercept of GLMM}

\item{b2}{Parameter related to between-subjects variance in the GLMM model}

\item{b3}{Parameter related to extra-dispersion in the GLMM model}

\item{b4}{Parameter related to extra zero probability in the GLMM model}

\item{sb}{Between-methods variability}
}
\value{
Scalar
}
\description{
ICC for Zero-Inflated Negative Binomial with proportional extradispersion. Concordance setting
}
\keyword{internal}
