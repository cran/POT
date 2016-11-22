\name{specdens}
\alias{specdens}
\title{Spectral Density Plot}
\description{
  Plot the spectral density for a bivariate extreme value distribution
  or an extreme Markov chain model. 
}
\usage{
specdens(fitted, main, plot = TRUE, \dots)
}
\arguments{
  \item{fitted}{An object of class \code{'bvpot'} or
    \code{'mcpot'}. Most often, the return object of the
    \code{\link{fitbvgpd}} or \code{\link{fitmcgpd}} function.}
  \item{main}{The title of the graphic window. May be missing.}
  \item{plot}{Logical. Should the spectral density be plotted? The
    default is to plot it.}
  \item{\dots}{Other options to be passed to the \code{\link{plot}}
    function.}
}
\value{
  Plot the spectral density for a fitted bivariate extreme value
  distribution. Moreover, the spectral density is returned invisibly.
}
\details{
  Any bivariate extreme value distribution has the following
  representation:

  \deqn{G(y_1, y_2) = \exp\left[ - \int_0^1 \max\left( \frac{q}{z_1},
    \frac{1-q}{z_2} \right) dH(q) \right]}{G(z_1, z_2) = exp{ -int_0^1
      max[q/z_1, (1-q)/z_2] dH(q)}}

  where \eqn{H}{H} holds:
  \deqn{\int_0^1 q dH(q) = \int_0^1 (1-q) dH(q) = 1}{int_0^1 q dH(q) =
    int_0^1 (1-q) dH(q) = 1}

  \eqn{H}{H} is called the spectral measure with density
  \eqn{h}{h}. Thus, \eqn{h}{h} is called the spectral density. In
  addition, \eqn{H}{H} has a total mass of 2.

  For two independent random variables, the spectral measure consists of
  two points of mass 1 at \eqn{q =0,1}{q=0,1}. For two perfect dependent
  random variables, the spectral measure consists of a single point of
  mass 2 at \eqn{q=0.5}{q=0.5}. 
}
\seealso{
  \code{\link{retlev.bvpot}}, \code{\link{pickdep}} and
  \code{\link{plot.bvpot}}
}
\examples{
par(mfrow=c(1,2))
##Spectral density for a Markov Model
mc <- simmc(1000, alpha = 0.25, model = "log")
mc <- qgpd(mc, 0, 1, 0.1)
Mclog <- fitmcgpd(mc, 0, "log")
specdens(Mclog)
##Spectral density for a bivariate POT model
x <- rgpd(500, 5, 1, -0.1)
y <- rgpd(500, 2, 0.2, -0.25)
Manlog <- fitbvgpd(cbind(x,y), c(5,2), "anlog")
specdens(Manlog)
}
\author{Mathieu Ribatet}
\keyword{hplot}
\keyword{models}
