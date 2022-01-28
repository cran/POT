\name{logLik.pot}
\alias{logLik.pot}
\title{Extract Log-Likelihood}
\description{Extract Log-Likelihood for object of class \code{'pot'}}
\usage{
\method{logLik}{pot}(object, \dots)
}
\arguments{
\item{object}{An object of class \code{'pot'}. Most often, this is an
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
fmle <- fitgpd(x, 0)
logLik(fmle)
}
\author{Mathieu Ribatet}
\keyword{models}
