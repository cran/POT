\name{print.bvpot}
\alias{print.bvpot}
\title{Printing bvpot objects}
\description{Print a ``bvpot'' object}
\usage{
\method{print}{bvpot}(x, digits = max(3, getOption("digits") - 3), \dots)
}
\arguments{
  \item{x}{An object of class ``bvpot''. Most often, returns of the
    \code{\link{fitbvgpd}} function.}
  \item{digits}{The number of digits to be printed.}
  \item{\dots}{Other options to be passed to the \code{\link{print}}
    function.}
}
\value{
  Print on screen.
}
\seealso{
  \code{\link{print.uvpot}}, \code{\link{print.mcpot}},
  \code{\link{print}} 
}
\examples{
set.seed(123)
x <- rgpd(500, 0, 1, 0.2)
y <- rgpd(500, 2, 0.5, -0.1)
Mlog <- fitbvgpd(cbind(x, y), c(0, 2))
Mlog
}
\author{Mathieu Ribatet}
\keyword{print}
