\name{fitpp}
\alias{fitpp}

\title{Fitting the point process characterisation to exceedances above a
  threshold}

\description{This function estimates the point process characterisation
  from exceedances above a threshold.}

\usage{fitpp(data, threshold, noy = length(data) / 365.25, start, \dots,
std.err.type = "observed", corr = FALSE, method = "BFGS", warn.inf = TRUE)}
\arguments{
  \item{data}{A numeric vector.}
  \item{threshold}{A numeric value giving the threshold for the
    GPD.}
  \item{noy}{Numeric. The number of year of observation.}
  \item{start}{A named list that gives the starting values for the
    optimization routine. Each list argument must correspond to one
    parameter to be estimated. May be missing.}
  \item{\dots}{Other optional arguments to be passed to the
    \code{\link{optim}} function, allow hand fixed parameters (only
    - see the Note section.}
  \item{std.err.type}{A character string. If "observed", the standard
    errors are derived from the observed Fisher information matrix. If
    "none", standard errors are not computed.}
  \item{corr}{Logical. Does the asymptotic correlation matrix has to be
    computed? Default is "not computed" - e.g. \code{FALSE}.}
  \item{method}{A character string specifying which numerical
    optimization procedure has to be used. See \code{\link{optim}} for
    more details.}
  \item{warn.inf}{Logical. If \code{TRUE} (default), users will be
    warned if the log-likelihood is not finite at starting values - as
    it may cause some problem during the optimation stage.}  
}
\value{
  This function returns a list with  components:
  \item{fitted.values}{A vector containing the estimated parameters.}
  \item{std.err}{A vector containing the standard errors.}
  \item{fixed}{A vector containing the parameters of the model that
    have been held fixed.}
  \item{param}{A vector containing all parameters (optimized and fixed).}
  \item{deviance}{The deviance at the maximum likelihood estimates.}
  \item{corr}{The correlation matrix.}
  \item{convergence, counts, message}{Components taken from the
    list returned by \code{\link{optim}} - for the \code{mle} method.}
  \item{threshold}{The threshold passed to argument \code{threshold}.}
  \item{nat, pat}{The number and proportion of exceedances.}
  \item{data}{The data passed to the argument \code{data}.}
  \item{exceed}{The exceedances, or the maxima of the clusters of
    exceedances.}
  \item{scale}{The scale parameter for the fitted generalized Pareto
    distribution.}
  \item{std.err.type}{The standard error type - for \code{'mle'}
    only. That is \code{Observed} Information matrix of Fisher.}
  \item{var.thresh}{Logical. Specify if the threshold is a varying one -
    \code{'mle'} only. For other methods, threshold is always constant
    i.e. \code{var.thresh = FALSE}. Not implemented yet.}
}
\references{
  Coles, S. (2001) \emph{An Introduction to Statistical
    Modelling of Extreme Values}. Springer Series in Statistics. London.

  Embrechts, P and Kluppelberg, C. and Mikosch, T (1997) \emph{Modelling
    Extremal Events for Insurance and Finance.} Springers.

  Pickands, J. (1975) Statistical Inference Using Extreme Order
  Statistics. Annals of Statistics. \bold{3}:119--131.
}

\examples{
x <- rgpd(1000, 0, 1, 0.2)
fitpp(x, 0)
}
\author{Mathieu Ribatet}
\keyword{htest}
