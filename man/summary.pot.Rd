\name{summary.pot}
\alias{summary.pot}
\title{Compactly display the structure}
\description{Compactly display the structure of an object of class \code{'pot'}}
\usage{
\method{summary}{pot}(object, \dots)
}
\arguments{
\item{object}{An object of class \code{'pot'}. Most often, this is an
    object return by the \code{\link{fitgpd}}, \code{\link{fitbvgpd}}
    and \code{\link{fitmcgpd}} functions.}
  \item{\dots}{Other arguments to be passed to the \code{\link{str}}
    function.}
}
\value{
  Standard \code{summary} object: see \code{\link{summary}}.
}
\seealso{
 \code{\link{summary}} 
}
\examples{
set.seed(123)
x <- rgpd(500, 0, 1, -0.15)
fmle <- fitgpd(x, 0)
summary(fmle)
}
\author{Christophe Dutang}
\keyword{models}
