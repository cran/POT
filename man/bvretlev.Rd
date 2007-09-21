\name{retlev.bvpot}
\alias{retlev.bvpot}
\title{Return Level Plot: Bivariate Case}
\description{
  Plot return levels for a fitted bivariate extreme value distribution.
}
\usage{
\method{retlev}{bvpot}(fitted, p = seq(0.75,0.95,0.05), main, n = 5000,
only.excess = FALSE, \dots)
}
\arguments{
  \item{fitted}{An object of class \code{"bvpot"}. Most often, the
    return object of the \code{\link{fitbvgpd}} function.}
  \item{p}{A vector of probabilities for which return levels must be
    drawn.}
  \item{main}{The title of the graphic window. May be missing.}
  \item{n}{The number (default: 5000) of points needed to draw return
    levels lines.}
  \item{only.excess}{Logical. If \code{FALSE} (the default), all
    observations are plotted, otherwise, only exceedances above at least
    one of the two thresholds are plotted.}
  \item{\dots}{Other parameters to pass to the \code{\link{plot}}
    function.}
}
\value{
  Plot return levels for a fitted bivariate extreme value
  distribution. Moreover, an invisible list is return which gives the
  points used to draw the current plot.
}
\details{
  Any bivariate extreme value distribution has the Pickands'
  representation form i.e.:
  \deqn{G(y_1, y_2) = \exp\left[ - \left(\frac{1}{z_1} + \frac{1}{z_2}
    \right) A( w ) \right]}{G(y1, y2) = exp( - (1/z1 + 1/z2) A(w) )}
  where \eqn{z_i}{zi} corresponds to \eqn{y_i}{yi} transformed to be
  unit Frechet distributed and \eqn{w = \frac{z_2}{z_1 + z_2}}{w = z2 /
    (z1 + z2)} which lies in \eqn{[0,1]}{[0,1]}.
  
  Thus, for a fixed probability \eqn{p}{p} and \eqn{w}{w}, we have the
  corresponding \eqn{z_1}{z1}, \eqn{z_2}{z2} values:
  \deqn{z_1 = - \frac{A(w)}{w \log(p)}}{z1 = - A(w) / (w log(p))}
  \deqn{z_2 = \frac{z_1 w}{1 - w}}{z2 = z1 w / (1-w)}
  
  At last, the \eqn{z_i}{zi} are transformed back to their original
  scale.  
}
\seealso{
  \code{\link{fitbvgpd}}, \code{\link{plot}}
}
\examples{
x <- rbvgpd(1000, alpha = 0.25, mar1 = c(0, 1, 0.25))
Mlog <- fitbvgpd(x, c(0, 0), "log")
retlev(Mlog)
}
\author{Mathieu Ribatet}
\keyword{hplot}
\keyword{models}
