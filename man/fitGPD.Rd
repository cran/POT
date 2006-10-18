\name{Fit the GP Distribution}
\alias{fitgpd}

\title{Fitting a GPD to Peaks Over a Threshold}
\description{Maximum Likelihood, Unbiased Probability Weigthed Moments,
  Biased Probability Weighted Moments and Moments Estimator to fit Peaks
  Over a Threshold to a GP distribution.}

\usage{fitgpd(data, threshold, method = "mle", \dots)}
\arguments{
  \item{data}{A numeric vector.}
  \item{threshold}{A numeric value giving the threshold for the
    GPD. \code{'mle'} method allows varying threshold; so that threshold
  could be for this case a numeric vector. Be careful, varying
  thresholds are used cyclically if length doesn't match with
  \code{data}.}
  \item{method}{A string giving the names of the estimator. It can be
    \code{'mle'} (the default), \code{'moments'}, \code{'pwmu'},
    \code{'pwmb'}, \code{'mdpd'}, \code{'med'} and \code{'pickands'} for
    the maximum likelihood, moments, unbiased probability weighted
    moments, biased probability weigthed moments, minimum density power
    divergence, medians and pickands' estimators respectively.}
  \item{\dots}{Other optional arguments to be passed to the
    \code{\link{optim}} function and allow hand fixed parameters.}
}
\value{
  This function returns a list with  composants:
  \item{fitted.values}{A vector containing the maximum likelihood estimates.}
  \item{std.err}{A vector containing the standard errors.}
  \item{fixed}{A vector containing the parameters of the model that
    have been held fixed.}
  \item{param}{A vector containing all parameters (optimized and fixed).}
  \item{deviance}{The deviance at the maximum likelihood estimates.}
  \item{corr}{The correlation matrix - for the \code{mle} method.}
  \item{convergence, counts, message}{Components taken from the
    list returned by \code{\link{optim}} - for the \code{mle} method.}
  \item{threshold}{The threshold passed to argument \code{threshold}.}
  \item{nat, pat}{The number and proportion of exceedences.}
  \item{data}{The data passed to the argument \code{data}.}
  \item{exceedances}{The exceedences, or the maxima of the clusters
    of exceedences.}
  \item{scale}{The scale parameter for the fitted generalized Pareto
    distribution.}
  \item{std.err.type}{The standard error type - for \code{'mle'}
    only. That is \code{Observed} or \code{Expected} Information matrix
    of Fisher.}
  \item{var.thresh}{Logical. Specify if the threshold is a varying one -
    \code{'mle'} only. For other methods, threshold is always constant
    i.e. \code{var.thresh = FALSE}.}
}
\note{The Maximum Likelihood estimator is obtained through a sligthly
  modified version of Alec Stephenson's \code{fpot.norm} function in the
  \code{evd} package.
}
\references{
  Coles, S. (2001) \emph{An Introduction to Statistical
    Modelling of Extreme Values}. Springer Series in Statistics. London.

  Hosking, J. and Wallis, J. (1987) Parameters and Quantile Estimation
  for the Generalized Pareto Distribution. Technometrics
  \bold{29}:339--349.

  Juarez, S. and Schucany, W. (2004) Robust and Efficient Estimation for
  the Generalized Pareto Distribution. Extremes \bold{7}:237--251.

  Peng, L. and Welsh, A. (2001) Robust Estimation of the Generalized
  Pareto Distribution. Extremes \bold{4}:53--65.

  Embrechts, P and Klüppelberg, C. and Mikosch, T (1997) \emph{Modelling
    Extremal Events for Insurance and Finance.} Springers.

  Pickands, J. (1975) Statistical Inference Using Extreme Order
  Statistics. Annals of Statistics. \bold{3}:119--131.  
}

\examples{
x <- rgpd(200, 1, 2, 0.25)
mle <- fitgpd(x, 1, "mle")$param
pwmu <- fitgpd(x, 1, "pwmu")$param
pwmb <- fitgpd(x, 1, "pwmb")$param
pickands <- fitgpd(x, 1, "pickands")$param    ##Check if Pickands estimates
                                              ##are valid or not !!!
med <- fitgpd(x, 1, "med", start = mle)$param ##Sometimes the fitting algo is not
                                              ##accurate. So specify
                                              ##good starting values is
                                              ##a good idea.  
mdpd <- fitgpd(x, 1, "mdpd")$param

print(rbind(mle, pwmu, pwmb, pickands, med, mdpd))
}
\author{Mathieu Ribatet}
\keyword{htest}
