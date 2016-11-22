\name{fitbvgpd}
\alias{fitbvgpd}
\title{Fitting Bivariate Peaks Over a Threshold Using Bivariate
  Extreme Value Distributions}
\description{
  Fitting a bivariate extreme value distribution to bivariate
  exceedances over thresholds using censored maximum likelihood
  procedure. 
}
\usage{
fitbvgpd(data, threshold, model = "log", start, \dots, cscale = FALSE,
cshape = FALSE, std.err.type = "observed", corr = FALSE, warn.inf = TRUE,
method = "BFGS")
}
\arguments{
  \item{data}{A matrix with two columns which gives the observation
    vector for margin 1 and 2 respectively. \code{NA} values are
    considered to fall below the threshold.}
  \item{threshold}{A numeric vector for the threshold (of length 2).}
  \item{model}{A character string which specifies the model used. Must
    be one of \code{log} (the default), \code{alog}, \code{nlog},
    \code{anlog}, \code{mix} and \code{amix} for the logistic, asymmetric
  logistic, negative logistic, asymmetric negative logistic, mixed and
  asymmetric mixed models.}
  \item{start}{Optional. A list for starting values in the fitting
    procedure.}
  \item{\dots}{Additional parameters to be passed to the
    \code{\link{optim}} function or to the bivariate model. In
    particular, parameter of the model can be hand fixed.}
  \item{cscale}{Logical. Should the two scale parameters be equal?}
  \item{cshape}{Logical. Should the two shape parameters be equal?}
  \item{std.err.type}{The type of the standard error. Currently, one
    must specify \code{"observed"} for observed Fisher information
    matrix or \code{"none"} for no computations of the standard
    errors.}
  \item{corr}{Logical. Should the correlation matrix be computed?}
  \item{warn.inf}{Logical. Should users be warned if likelihood is not
    finite at starting values?}
  \item{method}{The optimization method, see \code{\link{optim}}.}
}
\details{
  The bivariate exceedances are fitted using censored likelihood
  procedure. This methodology is fully described in Ledford (1996).

  Most of models are described in Kluppelberg (2006).
}
\section{Warnings}{
  Because of numerical problems, their exists artificial numerical
  constraints imposed on each model. These are:
  \itemize{
    \item For the logistic and asymmetric logistic models: \eqn{\alpha}{alpha}
    must lie in [0.05, 1] instead of [0,1];
    \item For the negative logistic model: \eqn{\alpha}{alpha} must lie in
    [0.01, 15] instead of \eqn{[0,\infty[}{[0, infinity[};
    \item For the asymmetric negative logistic model: \eqn{\alpha}{alpha} must
    lie in [0.2, 15] instead of \eqn{[0,\infty[}{[0, infinity[};
    \item For the mixed and asymmetric mixed models: None artificial numerical
    constraints are imposed.
  }

  For this purpose, users must check if estimates are near these
  artificial numerical constraints. Such cases may lead to substantial
  biases on the GP parameter estimates. One way to detect quickly if
  estimates are near the border constraints is to look at the standard
  errors for the dependence parameters. Small values (i.e. < 1e-5) often
  indicates that numerical constraints have been reached.

  In addition, users must be aware that the mixed and asymmetric mixed
  models can not deal with perfect dependence.

  Thus, user may want to plot the Pickands' dependence function to see
  if variable are near independence or dependence cases using the
  \code{\link{pickdep}} function.
}
\value{
  The function returns an object of class \code{c("bvpot","pot")}. As
  usual, one can extract several features using \code{\link{fitted}} (or
  \code{\link{fitted.values}}), \code{\link{deviance}},
  \code{\link{logLik}} and \code{\link{AIC}}  functions.  
  
  \item{fitted.values}{The maximum likelihood estimates of the bivariate
  extreme value distribution.}
  \item{std.err}{A vector containing the standard errors - only present
    when the observed information matrix is not singular.}
  \item{var.cov}{The asymptotic variance covariance matrix - only
    presents when the observed information matrix is not singular.}
  \item{deviance}{The deviance.}
  \item{corr}{The correlation matrix.}
  \item{convergence, counts, message}{Informations taken from the
    \code{\link{optim}} function.}
  \item{threshold}{The marginal thresholds.}
  \item{pat}{The marginal proportion above the threshold.}
  \item{nat}{The marginal number above the threshold.}
  \item{data}{The bivariate matrix of observations.}
  \item{exceed1, exceed2}{The marginal exceedances.}
  \item{call}{The call of the current function.}
  \item{model}{The model for the bivariate extreme value distribution.}
  \item{chi}{The chi statistic of Coles (1999). A value near 1 (resp. 0)
    indicates perfect dependence (resp. independence).}
}
\references{
  Coles, S., Heffernan, J. and Tawn, J. (1999) Dependence Measure for
  Extreme Value Analyses. Extremes, \bold{2}:4 339--365.
  
  Kl\"uppelberg, C., and May A. (2006) Bivariate extreme value
  distributions based on polynomial dependence functions. Mathematical
  Methods in the Applied Sciences, \bold{29}: 1467--1480.


  Ledford A., and Tawn, J. (1996) Statistics for near Independence in
  Multivariate Extreme Values. Biometrika, \bold{83}: 169--187.
}
\author{Mathieu Ribatet}
 
\seealso{
The following usual generic functions are available 
  \code{\link[POT:print.bvpot]{print}}, 
  \code{\link[POT:plot.bvpot]{plot}} and
  \code{\link[POT:anova.bvpot]{anova}}
  as well as new generic functions 
  \code{\link[POT:retlev.bvpot]{retlev}} and
  \code{\link[POT:convassess.bvpot]{convassess}}.
  
  See also \code{\link{pickdep}}, \code{\link{specdens}}.
  
  For optimization in R, see \code{\link{optim}}.
  
}


\examples{
x <- rgpd(1000, 0, 1, 0.25)
y <- rgpd(1000, 3, 1, -0.25)
ind <- fitbvgpd(cbind(x, y), c(0, 3), "log")
ind
ind2 <- fitbvgpd(cbind(x, y), c(0, 3), "log", alpha = 1)
ind2
ind3 <- fitbvgpd(cbind(x, y), c(0, 3), cscale = TRUE)
ind3
##The mixed model can not deal with perfect dependent variables
##Thus, there is a substantial bias in GPD parameter estimates
dep <- fitbvgpd(cbind(x, x + 3), c(0, 3), "mix")
dep
}
\keyword{models}
