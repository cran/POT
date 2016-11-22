\name{print.uvpot}
\alias{print.uvpot}
\title{Printing uvpot objects}
\description{Print an ``uvpot'' object}
\usage{
\method{print}{uvpot}(x, digits = max(3, getOption("digits") - 3), \dots)
}
\arguments{
\item{x}{An object of class \code{'uvpot'}. Most often, returns of the
    \code{\link{fitgpd}} function.}
  \item{digits}{The number of digits to be printed.}
  \item{\dots}{Other options to be passed to the \code{\link{print}}
    function.}
}
\value{
  Print on screen.
}
\seealso{
  \code{\link{print.bvpot}}, \code{\link{print.mcpot}},
  \code{\link{print}} 
}
\examples{
x <- rgpd(500, 0, 1, 0.2)
MLE <- fitgpd(x, 0)
MLE
}
\author{Mathieu Ribatet}
\keyword{print}
