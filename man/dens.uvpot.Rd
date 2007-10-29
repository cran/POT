\name{dens.uvpot}
\alias{dens.uvpot}
\title{Density Plot: Univariate Case}
\description{
  Density plot for univariate POT models. 
}
\usage{
\method{dens}{uvpot}(fitted, main, xlab, ylab, dens.adj = 1,
kern.lty = 2, rug = TRUE, plot.kernel = TRUE, plot.hist = TRUE,
hist.col = NULL, \dots)
}
\arguments{
  \item{fitted}{A object of class \code{``uvpot''}. Most often, the
    return of the \code{\link{fitgpd}} function.}
  \item{main}{The title of the graphic. If missing, the title is set to
    ``Density Plot''.}
  \item{xlab,ylab}{The labels for the x and y axis. If missing, they are
    set to ``Quantile'' and ``Density'' respectively.}
  \item{dens.adj}{Numeric. The adjustment for the kernel density
    estimation in the \code{\link{density}} function. The default is 1.}
  \item{kern.lty}{The line type for the kernel density estimation. This
    corresponds to the \code{``lty''} option of the \code{\link{lines}}
    functions. The default is 2.}
  \item{rug}{Logical. Should we call the \code{\link{rug}} function?
    Default is \code{TRUE}.}
  \item{plot.kernel}{Logical. Should the kernel density estimate be
    plotted?}
  \item{plot.hist}{Logical. Should the histogram be plotted?}
  \item{hist.col}{The color to fill the histogram.}
  \item{\dots}{Other arguments to be passed to the \code{\link{plot}}
    function.}
}
\value{
  A graphical window.
}
\details{
  The density plot consists of plotting on the same windows the
  theoretical density and a kernel estimation one. If the theoretical
  model is correct, then the two densities should be ``similar''.
}
\seealso{
  \code{\link{dens}}, \code{\link{dens.uvpot}}
}
\examples{
x <- rgpd(75, 1, 2, 0.1)
pwmu <- fitgpd(x, 1, "pwmu")
dens(pwmu)
}
\author{Mathieu Ribatet}
\keyword{hplot}
