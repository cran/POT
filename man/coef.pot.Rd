\name{coef.pot}
\alias{coef.pot}
\title{Extract model coefficients of a \code{'pot'} model}
\description{
\code{coef} extracts model coefficients of an object of class \code{'pot'}}
\usage{
\method{coef}{pot}(object, \dots)
}
\arguments{
\item{object}{An object of class \code{'pot'}. Most often, this is an
    object return by the \code{\link{fitgpd}}, \code{\link{fitbvgpd}}
    and \code{\link{fitmcgpd}} functions.}
  \item{\dots}{Other arguments to be passed to the \code{\link{str}}
    function.}
}
\value{
  Standard \code{coef} object: see \code{\link{coef}}.
}
\seealso{
 \code{\link{coef}} 
}
\examples{
set.seed(123)
x <- rgpd(500, 0, 1, -0.15)
fmle <- fitgpd(x, 0)
coef(fmle)
}
\author{Christophe Dutang}
\keyword{models}
