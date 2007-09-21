\name{logLik.pot}
\alias{logLik.pot}
\title{Extract Log-Likelihood}
\description{Extract Log-Likelihood for object of class ``pot''}
\usage{
\method{logLik}{pot}(object, \dots)
}
\arguments{
  \item{object}{An object of class ``pot''. Most often, this is an
    object return by the \code{\link{fitgpd}}, \code{\link{fitbvgpd}}
    and \code{\link{fitmcgpd}} functions.}
  \item{\dots}{Other arguments to be passed to the \code{\link{logLik}}
    function.}
}
\value{
  Standard \code{logLik} object: see \code{\link{logLik}}.
}
\seealso{
 \code{\link{logLik}} 
}
\examples{
x <- rgpd(500, 0, 1, -0.15)
mle <- fitgpd(x, 0)
logLik(mle)
}
\author{Mathieu Ribatet}
\keyword{models}
