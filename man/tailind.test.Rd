\name{tailind.test}
\alias{tailind.test}
\title{Testing for Tail Independence in Extreme Value Models}
\description{Several tests for tail independence (e.g. asymptotic
  independence) for a bivariate extreme value distribution
}
\usage{
tailind.test(data, c = -0.1, emp.trans = TRUE, chisq.n.class = 4)
}
\arguments{
  \item{data}{A matrix with two columns given the data.}
  \item{c}{A negative numeric. Must be close to zero to approximate
    accurately asymptotic results.}
  \item{emp.trans}{Logical. If \code{TRUE} (the default),
    \code{"data"} is transformed to reverse exponential using
    empirical estimates. Otherwise, \code{"data"} is supposed to be
    reverse exponential distributed.}
  \item{chisq.n.class}{A numeric given the number of classes for the Chi
  squared test.}
}
\value{
  This function returns a table with the Neymann-Pearson, Fisher,
  Kolmogorov-Smirnov and Chi-Square statistics and the related
  p-values.
}
\details{
  These tests are based on an asymptotic results shown by Falk and Michel
  (2006). Let \eqn{(X,Y)}{(X,Y)} be a random vector which follows in its
  upper tail a bivariate extreme value distribution with reverse
  exponential margins. The conditional distribution function of
  \eqn{X+Y}{X+Y}, given that \eqn{X+Y>c}{X+Y>c}, converges to
  \eqn{F(t)=t^2}{F(t)=t^2}, \eqn{t \in[0,1]}{t in [0,1]}, if \eqn{c
    \rightarrow 0^{-}}{t tends to 0} iff \eqn{X}{X} and \eqn{Y}{Y} are
  asymptotically independent. Otherwise, the limit is \eqn{F(t) =
    t}{F(t) = t}
}
\references{
  Falk, M. and Michel, Rene(2006) Testing for tail independence in
  extreme value models. \emph{Annals of the Institute of Statistical
    Mathematics} \bold{58}: 261--290
}
\seealso{
  \code{\link{chimeas}}, \code{\link{specdens}}
}
\examples{
##A total independence example
x <- rbvgpd(7000, alpha = 1, mar1 = c(0, 1, 0.25))
tailind.test(x)
##An asymptotically dependent example
y <- rbvgpd(7000, alpha = 0.75, model = "nlog", mar1 = c(0, 1, 0.25),
mar2 = c(2, 0.5, -0.15))
tailind.test(y)
##A perfect dependence example
z <- rnorm(7000)
tailind.test(cbind(z, 2*z - 5))
}
\author{Mathieu Ribatet}
\keyword{htest}