\name{plot.bvpot}
\alias{plot.bvpot}
\title{Graphical Diagnostics: the Bivariate Extreme Value
  Distribution Model.}
\description{
  Plot several graphics to judge goodness of fit of the fitted model.
}
\usage{
\method{plot}{bvpot}(x, mains, which = 1:3, ask = nb.fig < length(which)
&& dev.interactive(), \dots)
}
\arguments{
  \item{x}{An object of class \code{"bvpot"}. Most often, the object
    returned by the \code{\link{fitbvgpd}} function.}
  \item{mains}{May be missing. If present a 3--vector of character
    strings which gives the titles of the plots.}
  \item{which}{a numeric vector which specifies which plot must be drawn
    : \code{'1'} for Pickands' Dependence Function plot, \code{'2'} for
    a bivariate return level plot, \code{'3'} for the spectral density
    plot.}
  \item{ask}{Logical. If \code{TRUE}, user is asked before each plot.}
  \item{\dots}{Other parameters to pass to the \code{\link{plot}}
    function.}
}
\value{
  Several plots.
}
\seealso{
  \code{\link{fitbvgpd}}
}
\examples{
x <- rgpd(1000, 0, 1, 0.25)
y <- 3 + x + rnorm(1000)
Mlog <- fitbvgpd(cbind(x, y), c(0, 1), "log")
layout(matrix(c(1,1,2,2,0,3,3,0), 2, byrow = TRUE))
plot(Mlog)
}
\author{Mathieu Ribatet}
\keyword{hplot}
