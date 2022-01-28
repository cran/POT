\name{qq}
\alias{qq}
\alias{qq.uvpot}
\title{Quantile Quantile Plot}
\description{
\code{qq} is a generic function used to show quantile-quantile plot. 
The function invokes particular \code{methods} 
which depend on the \code{\link{class}} of the first argument.
So the function makes a quantile quantile plot for univariate POT models. 

}
\usage{
qq(object, \dots)

\method{qq}{uvpot}(object, main, xlab, ylab, ci = TRUE, \dots)
}
\arguments{
  \item{object}{A fitted object. When using the POT package, an object 
  of class \code{'uvpot'}. Most often, the
    return of the \code{\link{fitgpd}} function.}
  \item{main}{The title of the graphic. If missing, the title is set to
\code{"QQ-plot"}.}
  \item{xlab,ylab}{The labels for the x and y axis. If missing, they are
set to \code{"Model"} and \code{"Empirical"} respectively.}
  \item{ci}{Logical. If \code{TRUE} (the default), 95\% intervals are
    plotted.}
  \item{\dots}{Other arguments to be passed to the \code{\link{plot}}
    function.}
}
\value{
  A graphical window.
}
\details{
  The quantile quantile plot consists of plotting the observed quantiles
  in function of the theoretical ones. The theoretical quantiles
  \eqn{Q_{Theo, j}}{Q_{Theo, j}} are computed from the fitted GPD, that
  is:

  \deqn{Q_{Theo, j} = F^{-1}(p_j)}{Q_{Theo, j} = F^{-1}(p_j)}
  where \eqn{F^{-1}}{F^{-1}} is the fitted quantile function and
  \eqn{p_j}{p_j} are empirical probabilities defined by :

  \deqn{p_{j:n} = \frac{j - 0.35}{n}}{p_{j:n} = (j - 0.35) / n}
  where \eqn{n}{n} is the total number of observations - see Hosking
  (1995).

  If the theoretical model is correct, then points should be ``near''
  the line \eqn{y=x}{y=x}.
}
\seealso{
  \code{\link{qq}}, \code{\link{qq.uvpot}}
}
\examples{
x <- rgpd(75, 1, 2, 0.1)
pwmu <- fitgpd(x, 1, "pwmu")
qq(pwmu)
}
\references{
  Hosking, J. R. M. and Wallis, J. R. (1995). A comparison of unbiased
  and plotting-position estimators of L moments. Water Resources
  Research. \bold{31}(8): 2019--2025.
}
\author{Mathieu Ribatet}
\keyword{hplot}
