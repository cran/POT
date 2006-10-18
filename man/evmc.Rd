\name{simmc}
\alias{simmc}
\title{Simulate Markov Chains With Extreme Value
  Dependence Structures}
\description{
  Simulation of first order Markov chains, such that each pair
  of consecutive values has the dependence structure of one of
  nine parametric bivariate extreme value distributions. 
}
\usage{
simmc(n, alpha, model = "log", asCoef, asCoef1, asCoef2, margins =
"uniform")
}
\arguments{
  \item{n}{Number of observations.}
  \item{alpha}{Dependence parameter for the logistic, asymmetric
    logistic, negative logistic, asymmetric negative logistic, mixed and
    asymetric mixed models.}
  \item{asCoef,asCoef1,asCoef2}{The asymetric coefficients for the
    asymetric logisitc, asymetric negative logistique and asymetric
    mixed models.}
  \item{model}{The specified model; a character string. Must be
    either \code{"log"} (the default), \code{"alog"}, \code{"nlog"},
    \code{"anlog"}, \code{"mix"} or \code{"amix"}, for the logistic, asymmetric logistic,
    negative logistic, asymmetric negative logistic, mixed and asymmetric mixed
    models respectively.}
  \item{margins}{The marginal distribution of each value; a
    character string. Must be either \code{"uniform"} (the
    default), \code{"rweibull"}, \code{"frechet"} or
    \code{"gumbel"}, for the uniform, standard reversed Weibull,
    standard Gumbel and standard Frechet distributions respectively.}
}
\value{
  A numeric vector of length \code{n}.
}

\examples{
simmc(100, alpha = 0.1, model = "log")
simmc(100, alpha = 1.2, model = "nlog", margins = "gum")
}
\keyword{distribution}

\author{Alec Stephenson (modified for the POT package by Mathieu
  Ribatet)}


