\name{retlev.uvpot}
\alias{retlev.uvpot}
\title{Return Level Plot: Univariate Case}
\description{
  Return level plot for univariate POT models. 
}
\usage{
\method{retlev}{uvpot}(fitted, npy, main, xlab, ylab, xlimsup,
ci = TRUE, points = TRUE, \dots)
}
\arguments{
  \item{fitted}{A object of class \code{``uvpot''}. Most often, the
    return of the \code{\link{fitgpd}} function.}
  \item{npy}{The mean \bold{N}umber of events \bold{P}er \bold{Y}ear (or
    more generally per block).if missing, setting it to 1.}
  \item{main}{The title of the graphic. If missing, the title is set to
    ``Return Level Plot''.}
  \item{xlab,ylab}{The labels for the x and y axis. If missing, they are
    set to ``Return Period (Years)'' and ``Return Level'' respectively.}
  \item{xlimsup}{Numeric. The right limit for the x-axis. If missing, a
    suited value is computed.}
  \item{ci}{Logical. Should the 95\% pointwise confidence interval be
    plotted?}
  \item{points}{Logical. Should observations be plotted?}
  \item{\dots}{Other arguments to be passed to the \code{\link{plot}}
    function.}
}
\value{
  A graphical window. In addition, it returns invisibly the return level
  function.
}
\details{
  The return level plot consists of plotting the theoretical quantiles
  in function of the return period (with a logarithmic scale for the
  x-axis). For the definition of the return period see the
  \code{\link{prob2rp}} function. Thus, the return level plot consists
  of plotting the points defined by:

  \deqn{(T(p), F^{-1}(p))}{(T(p), F^{-1}(p))}
  where \eqn{T(p)}{T(p)} is the return period related to the non
  exceedance probability \eqn{p}{p}, \eqn{F^{-1}}{F^{-1}} is the
  fitted quantile function.

  If \code{points = TRUE}, the probabilities \eqn{p_j}{p_j} related to
  each observation are computed using the following plotting position
  estimator proposed by Hosking (1995):

  \deqn{p_j = \frac{j - 0.35}{n}}{p_j = (j - 0.35) / n}
  where \eqn{n}{n} is the total number of observations.

  If the theoretical model is correct, the points should be ``close''
  to the ``return level'' curve.
}
\seealso{
  \code{\link{retlev}}, \code{\link{retlev.uvpot}}
}
\examples{
x <- rgpd(75, 1, 2, 0.1)
pwmu <- fitgpd(x, 1, "pwmu")
rl.fun <- retlev(pwmu)
rl.fun(100)
}
\references{
  Hosking, J. R. M. and Wallis, J. R. (1995). A comparison of unbiased
  and plotting-position estimators of L moments. Water Resources
  Research. \bold{31}(8): 2019--2025.
}
\author{Mathieu Ribatet}
\keyword{hplot}
