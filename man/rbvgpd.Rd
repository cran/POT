\name{bvgpd}
\alias{rbvgpd}
\alias{pbvgpd}
\title{Parametric Bivariate GPD}
\description{
  Density, distribution function and random generation for six different
  parametric bivariate GPD
}
\usage{
rbvgpd(n, alpha, model = "log", asCoef, asCoef1, asCoef2, mar1 =
c(0,1,0), mar2 = mar1)
pbvgpd(q, alpha, model = "log", asCoef, asCoef1, asCoef2, mar1 =
c(0,1,0), mar2 = mar1, lower.tail = TRUE)
}
\arguments{
  \item{n}{The number of observations to be simulated.}
  \item{q}{A matrix or vector with two columns at which the distribution
    is computed.}
  \item{alpha}{Dependence parameter for the logistic, asymmetric
    logistic, negative logistic, asymmetric negative logistic, mixed and
    asymmetric mixed models.}
  \item{model}{The specified model; a character string. Must be
    either \code{"log"} (the default), \code{"alog"}, \code{"nlog"},
    \code{"anlog"}, \code{"mix"} or \code{"amix"}, for the logistic,
    asymmetric logistic, negative logistic, asymmetric negative
    logistic, mixed and asymmetric mixed models respectively.}
  \item{asCoef,asCoef1,asCoef2}{The asymmetric coefficients for the
    asymmetric logistic, asymmetric negative logistic and asymmetric
    mixed models.}
  \item{mar1,mar2}{Vectors of length 3 giving the marginal parameters.}
  \item{lower.tail}{Logical. If \code{TRUE} (the default), P[X
    \eqn{\leq}{<=} x] is computed, else P[X \eqn{\geq}{>=} x].}
}
\value{
  Generate a random vector of length \code{n}.
}
\details{
  The logistic and asymmetric logistic models respectively are
  simulated using bivariate versions of Algorithms 1.1 and 1.2 in
  Stephenson(2003).
  All other models are simulated using a root finding algorithm
  to simulate from the conditional distributions.
}
\references{
  Stephenson, A. G. (2003)
  Simulating multivariate extreme value distributions of logistic type.
  \emph{Extremes}, \bold{6}(1), 49--60.
}
\examples{
x <- rbvgpd(1000, alpha = 0.25, model = "log", mar1 = c(0,1,0.25), mar2
= c(2,0.5, -0.15))
y <- rbvgpd(1000, alpha = 0.75, model = "nlog", mar1 = c(0,1,0.25), mar2
= c(2,0.5, -0.15))
par(mfrow=c(1,2))
plot(x);plot(y)
}
\author{Mathieu Ribatet (Alec Stephenson for the C codes)}
\keyword{distribution}
