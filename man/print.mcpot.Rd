\name{print.mcpot}
\alias{print.mcpot}
\title{Printing mcpot objects}
\description{Print an ``mcpot'' object}
\usage{
\method{print}{mcpot}(x, digits = max(3, getOption("digits") - 3), \dots)
}
\arguments{
\item{x}{An object of class \code{'mcpot'}. Most often, returns of the
    \code{\link{fitmcgpd}} function.}
  \item{digits}{The number of digits to be printed.}
  \item{\dots}{Other options to be passed to the \code{\link{print}}
    function.}
}
\value{
  Print on screen.
}
\seealso{
  \code{\link{print.uvpot}}, \code{\link{print.bvpot}},
  \code{\link{print}} 
}
\examples{
x <- simmc(1000, alpha = 0.5)
x <- qgpd(x, 0, 1, 0.15)
Mc <- fitmcgpd(x, 0)
Mc
}
\author{Mathieu Ribatet}
\keyword{print}
