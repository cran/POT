\name{Fit the GP Distribution}
\alias{fitgpd}

\title{Fitting a GPD to Peaks Over a Threshold}

\description{Maximum (Penalized) Likelihood, Unbiased Probability
  Weighted Moments,Biased Probability Weighted Moments, Moments,
  Pickands', Minimum Density Power Divergence, Medians, Likelihood
  Moment and Maximum Goodness-of-Fit Estimators to fit Peaks Over a
  Threshold to a GP distribution.}

\usage{fitgpd(data, threshold, est = "mle", \dots)}
\arguments{
  \item{data}{A numeric vector.}
  \item{threshold}{A numeric value giving the threshold for the
    GPD. The \code{'mle'} estimator allows varying threshold; so that threshold
  could be for this case a numeric vector. Be careful, varying
  thresholds are used cyclically if length doesn't match with
  \code{data}.}
  \item{est}{A string giving the names of the estimator. It can be
    \code{'mle'} (the default), \code{'mple'}, \code{'moments'}, \code{'pwmu'},
    \code{'pwmb'}, \code{'mdpd'}, \code{'med'}, \code{'pickands'},
    \code{'lme'} and \code{'mgf'} for the maximum likelihood, maximum
    penalized likelihood, moments, unbiased
    probability weighted moments, biased probability weigthed moments,
    minimum density power divergence, medians, Pickands', likelihood
    moment and maximum goodness-of-fit estimators respectively.} 
  \item{\dots}{Other optional arguments to be passed to the
    \code{\link{optim}} function, allow hand fixed parameters (only
    for the \code{mle}, \code{mple} and \code{mgf} estimators) or passed
    several options to specific estimators - see the Note section.}
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
    only. That is \code{Observed} or \code{Expected} Information matrix
    of Fisher.}
  \item{var.thresh}{Logical. Specify if the threshold is a varying one -
    \code{'mle'} only. For other methods, threshold is always constant
    i.e. \code{var.thresh = FALSE}.}
}
\note{The Maximum Likelihood estimator is obtained through a slightly
  modified version of Alec Stephenson's \code{fpot.norm} function in the
  \code{evd} package.

  For the \code{'mple'} estimator, the likelihood function is penalized
  using the following function :

  \deqn{P(\xi) = \left\{
    \begin{array}{ll}
    1, & \xi \leq 0\\
    \exp\left[-\lambda \left(\frac{1}{1-\xi} - 1\right)^\alpha \right], &
    0 < \xi <1\\
    0, & \xi \geq 1
    \end{array}
  \right.}{P(xi)=1 if xi \le 0, exp[-lambda(1/(1-xi) -
    1)^alpha], if 0 < xi < 1, 0 if xi \ge 1}
  where \eqn{\alpha}{alpha} and \eqn{\lambda}{lambda} are the penalizing
  constants. Coles and Dixon (1999) suggest the use of
  \eqn{\alpha=\lambda=1}{alpha=lambda=1}.

  The \code{'lme'} estimator has a special parameter \code{'r'}. Zhang
  (2007) shows that a value of -0.5 should be accurate in most of the
  cases. However, other values such as \code{r < 0.5} can be
  explored. In particular, if \code{r} is approximatively equal to the
  opposite of the true shape parameter value, then the \code{lme}
  estimate is equivalent to the \code{mle} estimate.

  The \code{'pwmb'} estimator has special parameters \code{'a'} and
  \code{'b'}. These parameters are called the "plotting-position"
  values. Hosking and Wallis (1987) recommend the use of \code{a = 0.35}
  and \code{b = 0} (the default). However, different values can be
  tested.

  For the \code{'pwmu'} and \code{'pwmb'} approaches, one can pass the
  option \code{'hybrid = TRUE'} to use hybrid estimators as proposed by
  Dupuis and Tsao (1998). Hybrid estimators avoid to have no feasible
  points.

  The \code{mdpd} estimator has a special parameter \code{'a'}. This is
  a parameter of the "density power divergence". Juarez and Schucany
  (2004) recommend the use of \code{a = 0.1}, but any value of
  \code{a} such as \code{a > 0} can be used (small values are recommend
  yet).

  The \code{med} estimator admits two extra arguments \code{tol} and
  \code{maxit} to control the "stopping-rule" of the optimization
  process.

  The \code{mgf} approach uses goodness-of-fit statistics to estimate
  the GPD parameters. There are currently 8 different statitics: the
  Kolmogorov-Smirnov \code{"KS"}, Cramer von Mises \code{"CM"}, Anderson
  Darling \code{"AD"}, right tail Anderson Darling \code{"ADR"}, left
  tail Anderson Darling \code{"ADL"}, right tail Anderson Darling
  (second degree) \code{"AD2R"}, left tail Anderson Darling (second
  degree) \code{"AD2L"} and the Anderson Darling (second degree)
  \code{"AD2"} statistics.
}
\references{
  Coles, S. (2001) \emph{An Introduction to Statistical
    Modelling of Extreme Values}. Springer Series in Statistics. London.

  Coles, S. and Dixon, M. (1999) Likelihood-Based Inference for Extreme
  Value Models. Extremes \bold{2}(1):5--23.

  Dupuis, D. and Tsao (1998) M. A hybrid estimator for generalized
  Pareto and extreme-value distributions. Communications in
  Statistics-Theory and Methods \bold{27}:925--941.

  Hosking, J. and Wallis, J. (1987) Parameters and Quantile Estimation
  for the Generalized Pareto Distribution. Technometrics
  \bold{29}:339--349.

  Juarez, S. and Schucany, W. (2004) Robust and Efficient Estimation for
  the Generalized Pareto Distribution. Extremes \bold{7}:237--251.

  Luceno, A. (2006) Fitting the generalized Pareto distribution to data
  using maximum goodness-of-fit estimators. Computational Statistics and
  Data Analysis \bold{51}:904--917.

  Peng, L. and Welsh, A. (2001) Robust Estimation of the Generalized
  Pareto Distribution. Extremes \bold{4}:53--65.

  Embrechts, P and Kluppelberg, C. and Mikosch, T (1997) \emph{Modelling
    Extremal Events for Insurance and Finance.} Springers.

  Pickands, J. (1975) Statistical Inference Using Extreme Order
  Statistics. Annals of Statistics. \bold{3}:119--131.

  Zhang, J. (2007) Likelihood Moment Estimation for the Generalized
  Pareto Distribution. Australian and New Zealand Journal of
  Statistics. \bold{49}(1):69--77.
}

\examples{
x <- rgpd(200, 1, 2, 0.25)
mle <- fitgpd(x, 1, "mle")$param
pwmu <- fitgpd(x, 1, "pwmu")$param
pwmb <- fitgpd(x, 1, "pwmb")$param
pickands <- fitgpd(x, 1, "pickands")$param    ##Check if Pickands estimates
                                              ##are valid or not !!!
med <- fitgpd(x, 1, "med",                    ##Sometimes the fitting algo is not
start = list(scale = 2, shape = 0.25))$param  ##accurate. So specify
                                              ##good starting values is
                                              ##a good idea.  
mdpd <- fitgpd(x, 1, "mdpd")$param
lme <- fitgpd(x, 1, "lme")$param
mple <- fitgpd(x, 1, "mple")$param
ad2r <- fitgpd(x, 1, "mgf", stat = "AD2R")$param

print(rbind(mle, pwmu, pwmb, pickands, med, mdpd, lme,
 mple, ad2r))

##Use PWM hybrid estimators
fitgpd(x, 1, "pwmu", hybrid = FALSE)

##Now fix one of the GPD parameters
##Only the MLE, MPLE and MGF estimators are allowed !
fitgpd(x, 1, "mle", scale = 2)
fitgpd(x, 1, "mple", shape = 0.25) 
}
\author{Mathieu Ribatet}
\keyword{htest}
