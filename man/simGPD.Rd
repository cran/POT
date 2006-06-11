\name{Generalized Pareto}
\alias{rgpd}
\alias{dgpd}
\alias{qgpd}
\alias{pgpd}

\title{The Generalized Pareto Distribution}
\description{Density, distribution function, quantile function and random
  generation for the GP distribution with location equal to 'loc',
  scale equal to 'scale' and shape equal to 'shape'.}

\usage{
rgpd(n, loc = 0, scale = 1, shape = 0)
pgpd(q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
qgpd(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
dgpd(x, loc = 0, scale = 1, shape = 0, log = FALSE)
}

\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.}
  \item{loc}{vector of the location parameters.}
  \item{scale}{vector of the scale parameters.}
  \item{shape}{a numeric of the shape parameter.}
  \item{lower.tail}{logical; if TRUE (default), probabilities are \eqn{P[ X
      \le x]}, otherwise, \eqn{P[X > x]}.}
  \item{log}{logical; if TRUE, probabilities p are given as log(p).}
}
\value{
  If 'loc', 'scale' and 'shape' are not specified they assume the default
  values of '0', '1' and '0', respectively.
  
  The GP distribution function for loc = \eqn{\mu}, scale = \eqn{\sigma}
  and shape = \eqn{\xi} is
  
  \deqn{G(z) = 1 - \left[ 1 + \frac{\xi ( z - \mu )}{ \sigma } \right] ^ { - 1 /
      \xi}}{%
    G(z) = 1 - [ 1 + \xi ( z - \mu ) / \sigma ]^(-1/\xi)}
  
  for \eqn{1 + \xi ( z - \mu ) / \sigma > 0} and \eqn{z > \mu}, where
  \eqn{\sigma > 0}. If
  \eqn{\xi = 0}, the distribution is defined by continuity corresponding to
  the exponential distribution. 
}
\examples{
dgpd(0.1)
rgpd(100, 1, 2, 0.2)
qgpd(seq(0.1, 0.9, 0.1), 1, 0.5, -0.2)
pgpd(12.6, 2, 0.5, 0.1)
}
\keyword{distribution}
