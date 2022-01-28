\name{fitmcgpd}
\alias{fitmcgpd}
\title{Fitting Markov Chain Models to Peaks Over a Threshold}
\description{
  Fitting a Markov chain to cluster exceedances using a bivariate
  extreme value distribution and a censored maximum likelihood
  procedure.
}
\usage{
fitmcgpd(data, threshold, model = "log", start, \dots, std.err.type =
"observed", corr = FALSE, warn.inf = TRUE, method = "BFGS")
}
\arguments{
  \item{data}{A vector of observations.}
  \item{threshold}{The threshold value.}
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
  \item{std.err.type}{The type of the standard error. Currently, one
    must specify \code{``observed''} for observed Fisher information
    matrix or \code{``none''} for no computations of the standard
    errors.}
  \item{corr}{Logical. Should the correlation matrix be computed?}
  \item{warn.inf}{Logical. Should users be warned if likelihood is not
    finite at starting values?}
  \item{method}{The optimization method, see \code{\link{optim}}.}
}
\value{
  The function returns an object of class \code{c("mcpot", "uvpot",
    "pot")}.  As usual, one can extract several features using
  \code{\link{fitted}} (or \code{\link{fitted.values}}),
  \code{\link{deviance}}, \code{\link{logLik}} and \code{\link{AIC}}
  functions.

   \item{fitted.values}{The maximum likelihood estimates of the Markov
  chain including estimated parameters of the bivariate extreme value
  distribution.}
  \item{std.err}{A vector containing the standard errors - only present
    when the observed information matrix is not singular.}
  \item{var.cov}{The asymptotic variance covariance matrix - only
    presents when the observed information matrix is not singular.}
  \item{deviance}{The deviance.}
  \item{corr}{The correlation matrix.}
  \item{convergence, counts, message}{Informations taken from the
    \code{\link{optim}} function.}
  \item{threshold}{The threshold.}
  \item{pat}{The proportion above the threshold.}
  \item{nat}{The number above the threshold.}
  \item{data}{The observations.}
  \item{exceed}{The exceedances.}
  \item{call}{The call of the current function.}
  \item{model}{The model for the bivariate extreme value distribution.}
  \item{chi}{The chi statistic of Coles (1999). A value near 1 (resp. 0)
    indicates perfect dependence (resp. independence).}
}
\details{
  The Markov Chain model is defined as follows:
  \deqn{L\left(y;\theta_1,\theta_2\right) = f\left(x_1; \theta_1\right)
    \prod_{i=2}^n f\left(y_i |
    y_{i-1};\theta_1,\theta_2\right)}{f(x_1;theta_1) f(x_2 | x_1;theta_1,
    theta_2) \dots f(x_n | x_{n-1}; theta_1, theta_2)}
  
  As exceedances above a (high enough) threshold are of interest, it is
  assumed that the marginal are GPD distributed, while the joint
  distribution is represented by a bivariate extreme value
  distribution. Smith et al. (1997) present theoretical results about
  this Markov Chain model.
  
  The bivariate exceedances are fitted using censored likelihood
  procedure. This methodology is fully described in Ledford (1996).

  Most of models are described in Kluppelberg (2006).
}
\section{Warnings}{
  Because of numerical problems, there exists artificial numerical
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

  In addition, we recommend to fix the marginal parameters. Indeed, even
  this is a two steps optimization procedure, this avoid numerical
  troubles - the likelihood function for the Markov chain model seems to
  be problematic. Thus, estimates are often better using the two stages
  approach.
}
\references{
  Kl\"uppelberg, C., and May A. (2006) Bivariate extreme value
  distributions based on polynomial dependence functions. Mathematical
  Methods in the Applied Sciences, \bold{29} 1467--1480.


  Ledford A., and Tawn, J. (1996) Statistics for near Independence in
  Multivariate Extreme Values. Biometrika, \bold{83} 169--187.


  Smith, R., and Tawn, J., and Coles, S. (1997) Markov chain models for
  threshold exceedances. Biometrika, \bold{84} 249--268
}
\seealso{
The following usual generic functions are available 
  \code{\link[POT:print.mcpot]{print}}, 
  \code{\link[POT:plot.mcpot]{plot}}
  as well as new generic functions 
  \code{\link[POT:retlev.mcpot]{retlev}} and
  \code{\link[POT:convassess.mcpot]{convassess}}.
  
  See also \code{\link{pickdep}}.
  
  For optimization in R, see \code{\link{optim}}.
  
}

\examples{
mc <- simmc(1000, alpha = 0.25)
mc <- qgpd(mc, 0, 1, 0.25)
##A first application when marginal parameter are estimated
fitmcgpd(mc, 0)

##Another one where marginal parameters are fixed
fmle <- fitgpd(mc, 0)
fitmcgpd(mc, 0, scale = fmle$param["scale"], shape = fmle$param["shape"])
}
\author{Mathieu Ribatet}
\keyword{models}
