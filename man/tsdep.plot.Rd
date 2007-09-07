\name{tsdep.plot}
\alias{tsdep.plot}
\title{Diagnostic for Dependence within Time Series Extremes}
\description{
  A diagnostic tool to assess for short range asymptotic dependence
  within a stationary time series.
}
\usage{
tsdep.plot(data, u, \dots, xlab, ylab, n.boot = 100, show.lines = TRUE,
lag.max, ci = 0.95, block.size = 5 * lag.max, angle = 90, arrow.length =
0.1)
}
\arguments{
  \item{data}{The time series observations.}
  \item{u}{The threshold.}
  \item{\dots}{Optional arguments to be passed to the \code{\link{plot}}
    function.}
  \item{xlab,ylab}{The x and y-axis labels.}
  \item{n.boot}{Numeric. The number of replicates to compute the
    bootstrap confidence interval.}
  \item{show.lines}{Logical. If \code{TRUE} (the default), the
    theoretical lines for the asymptotic dependence and ``near''
    independence are drawn.}
  \item{lag.max}{The maximum lag to be explored - may be missing.}
  \item{ci}{The level for the bootstrap confidence interval. The default
    is the 95\% confidence interval.}
  \item{block.size}{The size for the contiguous bootstrap approach.}
  \item{angle}{The angle at the end of the error bar. If \code{0}, error
    bars are only segments.}
  \item{arrow.length}{The length to be passed in the function
    \code{\link{arrows}}.}  
}
\value{
  This function plot the \eqn{\Lambda_\tau}{Lambda_tau} statictics
  against the lag. Bootstrap confidence intervals are also drawn. The
  function returns invisibly this statistic and the confidence bounds.
}
\details{
  Let \code{X_t} be a stationary sequence of unit Frechet random
  variables. By stationarity, the joint survivor function
  \eqn{\overline{F}_\tau(\cdot, \cdot)}{Fbar_tau(.,.)} of \eqn{(X_t,
    X_{t+\tau})}{(X_t, X_{t+tau})} does not depend on \eqn{t}{t}.

  One parametric representation for \eqn{\overline{F}_\tau(\cdot,
    \cdot)}{Fbar_tau(.,.)} is given by
  \deqn{\overline{F}_\tau(s,s)=L_\tau(s) s^{-1/\eta_\tau}}{Fbar_tau(s,s)
    = L_tau(x) s^{-1/eta_tau}}
  for some parameter \eqn{\eta_\tau \in (0,1]}{eta_tau in (0,1]} and a
  slowly varying function \eqn{L_\tau}{L_tau}.

  The \eqn{\Lambda_\tau}{Lambda_tau} statistic is defined by
  \deqn{\Lambda_\tau = 2 \eta_\tau - 1}{Lambda_tau = 2 eta_tau - 1}
  This statistic belongs to (-1,1] and is a measure of extremal
  dependence. \eqn{\Lambda_\tau = 1}{Lambda_tau = 1} corresponds to
  asymptotic dependence, \eqn{0 < \Lambda_\tau < 1}{0 < Lambda_tau < 1}
  to positive extremal association, \eqn{\Lambda_\tau = 0}{Lambda_tau =
    0} to ``near'' independence and \eqn{\Lambda_\tau < 0}{Lambda_tau <
    0} to negative extremal association.
}
\references{
  Ledford, A. and Tawn, J. (2003) Diagnostics for dependence within time
  series extremes. \emph{L. R. Statist. Soc. B.} \bold{65}, Part 2,
  521--543.

  Ledford, A. and Tawn, J (1996) Statistics for near independence in
  multivariate extreme values. \emph{Biometrika} \bold{83} 169--187.
}
\seealso{
  \code{\link{chimeas}}, \code{\link{tailind.test}}
}
\examples{
##An independent case
tsdep.plot(runif(5000), u = 0.95, lag.max = 5)
##Asymptotic dependence
mc <- simmc(5000, alpha = 0.2)
tsdep.plot(mc, u = 0.95, lag.max = 5)
}
\author{Mathieu Ribatet}
\keyword{htest}
