\name{Threshold Selection}
\alias{tcplot}
\title{Threshold Choice Plot}
\description{
  Plots of parameter estimates at various thresholds for peaks
  over threshold modelling, using the Generalized Pareto or Point
  Process representation. 
}
\usage{
tcplot(data, tlim, cmax = FALSE, r = 1,
    ulow = -Inf, rlow = 1, nt = 25, which = 1:npar, conf = 0.95,
    lty = 1, lwd = 1, type = "b", cilty = 1, ask = nb.fig <
    length(which) && dev.interactive(), \dots)
}
\arguments{
  \item{data}{A numeric vector.}
  \item{tlim}{A numeric vector of length two, giving the limits for the
    thresholds at which the model is fitted.}
  \item{cmax}{Logical; if \code{FALSE} (the default), the models are
    fitted using all exceedences over the thresholds. If \code{TRUE},
    the models are fitted using cluster maxima.}
  \item{r, ulow, rlow}{Arguments used for the identification of clusters
    of exceedences. Ignored if \code{cmax} is \code{FALSE} (the
    default).}
  \item{nt}{The number of thresholds at which the model is fitted.}
  \item{which}{If a subset of the plots is required, specify a
    subset of the numbers \code{1:npar}, where \code{npar} is
    the number of parameters.}
  \item{conf}{The (pointwise) confidence coefficient for the
    plotted confidence intervals. Use zero to suppress.}
  \item{lty, lwd}{The line type and width of the line connecting
    the parameter estimates.}
  \item{type}{The form taken by the line connecting the parameter
    estimates and the points denoting these estimates. Possible
    values include \code{"b"} (the default) for points joined by
    lines, \code{"o"} for overplotted points and lines, and
    \code{"l"} for an unbroken line with no points.}
  \item{cilty}{The line type of the lines depicting the
    confidence intervals.}
  \item{ask}{Logical; if \code{TRUE}, the user is asked before
    each plot.}
  \item{\dots}{Other arguments to be passed to the model fit
    function \code{fitgpd}.}
}
\details{
  For each of the \code{nt} thresholds a peaks over threshold model
  is fitted using the function \code{fitgpd}. The maximum likelihood
  estimates for the shape and the modified scale parameter (modified by
  subtracting the shape multiplied by the threshold) are plotted against
  the thresholds. If the threshold \code{u} is a valid threshold to be
  used for peaks over threshold modelling, the parameter estimates
  depicted should be approximately constant above \code{u}.
}
\value{
  A list is invisibly returned. Each component is a matrix with
  three columns giving parameter estimates and confidence limits.
}
\author{Stuart Coles and Alec Stephenson}

\seealso{\code{\link{fitgpd}}, \code{\link{mrlplot}}}
\examples{
data(ardieres)
flows <- ardieres[,"flows"]
tcplot(flows, tlim = c(0,15) )
}
\references{
  Coles, S. (2001). \emph{An Introduction to Statistical
    Modelling of Extreme Values}. Springer Series in Statistics. London.
  }
\keyword{hplot}
