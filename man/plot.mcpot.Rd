\name{plot.mcpot}
\alias{plot.mcpot}
\title{Graphical Diagnostics: Markov Chains for All Exceedances.}
\description{
  Plot several graphics to judge goodness of fit of the fitted model.
}
\usage{
\method{plot}{mcpot}(x, opy, exi, mains, which = 1:4, ask = nb.fig <
length(which) && dev.interactive(), acf.type = "partial", \dots)
}
\arguments{
  \item{x}{An object of class \code{"bvpot"}. Most often, the object
    returned by the \code{\link{fitbvgpd}} function.}
  \item{opy}{Numeric. The number of \bold{O}bservation \bold{P}er
    \bold{Y}ear (or more generally per block). If missing, the function
    warns and set it to 365.}
  \item{exi}{Numeric. The extremal index value. If missing, the
    estimator of Ferro and Segers (2003) is used.}
  \item{mains}{May be missing. If present a 4--vector of character
    strings which gives the titles of the plots.}
  \item{which}{a numeric vector which specifies which plot must be
    drawn: \code{'1'} for the auto correlation plot, \code{'2'} for
    Pickands' Dependence Function plot, \code{'3'} for the spectral
    density plot and \code{'4'} for a bivariate return level plot.}
  \item{ask}{Logical. If \code{TRUE}, user is asked before each plot.}
  \item{acf.type}{The type of auto correlation to be plotted. Must be
    one of \code{"correlation"}, \code{"covariance"} or \code{"partial"}
    (the default). See the \code{\link{acf}} function.}
  \item{\dots}{Other parameters to pass to the \code{\link{plot}}
    function.}
}
\value{
  Several plots and returns invisibly the return level function.
}
\note{
  For the return level plot, the observations are not plotted as these
  are dependent realisations. In particular, the return periods computed
  using the \code{\link{prob2rp}} are inaccurate.  
}
\seealso{
  \code{\link{fitmcgpd}}, \code{\link{acf}}, \code{\link{retlev}}
}
\examples{
set.seed(123)
mc <- simmc(200, alpha = 0.5)
mc <- qgpd(mc, 0, 1, 0.25)
Mclog <- fitmcgpd(mc, 1)
par(mfrow=c(2,2))
rlMclog <- plot(Mclog)
rlMclog(T = 3)
}
\section{Warning}{
  See the warning for the return level estimation in documentation of
  the \code{\link{retlev.mcpot}} function.
}
\references{
  Ferro, C. and Segers, J. (2003). Inference for clusters of extreme
  values. Journal of the Royal Statistical Society B. \bold{65}:
  545--556.
}
\author{Mathieu Ribatet}
\keyword{hplot}
