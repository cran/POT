\name{retlev}
\alias{retlev}
\alias{retlev.mcpot}
\alias{retlev.uvpot}
\title{Return Level Plot}
\description{
\code{retlev} is a generic function used to show return level plot. 
The function invokes particular \code{methods} 
which depend on the \code{\link{class}} of the first argument.
So the function makes a return level plot for POT models. 
}
\usage{
retlev(object, \dots)

\method{retlev}{uvpot}(object, npy, main, xlab, ylab, xlimsup,
ci = TRUE, points = TRUE, \dots)
\method{retlev}{mcpot}(object, opy, exi, main, xlab, ylab, xlimsup,
\dots)
}
\arguments{
  \item{object}{A fitted object. When using the POT package, an object 
  of class \code{'uvpot'} or \code{'mcpot'}. Most often, the
    return of  \code{\link{fitgpd}} or \code{\link{fitmcgpd}} functions.}
  \item{npy}{The mean \bold{N}umber of events \bold{P}er \bold{Y}ear (or
    more generally per block).if missing, setting it to 1.}
  \item{main}{The title of the graphic. If missing, the title is set to
\code{"Return Level Plot"}.}
  \item{xlab,ylab}{The labels for the x and y axis. If missing, they are
set to \code{"Return Period (Years)"} and \code{"Return Level"} respectively.}
  \item{xlimsup}{Numeric. The right limit for the x-axis. If missing, a
    suited value is computed.}
  \item{ci}{Logical. Should the 95\% pointwise confidence interval be
    plotted?}
  \item{points}{Logical. Should observations be plotted?}
  \item{\dots}{Other arguments to be passed to the \code{\link{plot}}
    function.}
  \item{opy}{The number of \bold{O}bservations \bold{P}er \bold{Y}ear (or
    more generally per block). If missing, it is set it to 365 i.e. daily
    values with a warning.}
  \item{exi}{Numeric. The extremal index. If missing, an estimate is
    given using the \code{\link{fitexi}} function.}
}
\value{
  A graphical window. In addition, it returns invisibly the return level
  function.
}
\details{

For class \code{"uvpot"}, the return level plot consists of plotting the theoretical quantiles
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
  
For class \code{"mcpot"}, let \eqn{X_1, \ldots,X_n}{X_1, \ldots,X_n} be the first \eqn{n}{n}
  observations from a stationary sequence with marginal distribution
  function \eqn{F}{F}. Thus, we can use the following (asymptotic)
  approximation:

  \deqn{\Pr\left[\max\left\{X_1,\ldots,X_n\right\} \leq x \right] =
    \left[ F(x) \right]^{n \theta}}{Pr[max{X_1,\ldots,X_n} <= x] =
    [F(x)]^(n theta)}

  where \eqn{\theta}{theta} is the extremal index.

  Thus, to obtain the T-year return level, we equate this equation to
  \eqn{1 - 1/T}{1 - 1/T} and solve for \eqn{x}{x}.  
}
\section{Warning}{
  For class \code{"mcpot"}, though this is computationally expensive, we recommend to give the
  extremal index estimate using the \code{\link{dexi}} function. Indeed,
  there is a severe bias when using the Ferro and Segers (2003)
  estimator - as it is estimated using observation and not the Markov
  chain model.
}

\seealso{
  \code{\link{prob2rp}}, \code{\link{fitexi}}.
}
\examples{

#for uvpot class
x <- rgpd(75, 1, 2, 0.1)
pwmu <- fitgpd(x, 1, "pwmu")
rl.fun <- retlev(pwmu)
rl.fun(100)

#for mcpot class
data(ardieres)
Mcalog <- fitmcgpd(ardieres[,"obs"], 5, "alog")
retlev(Mcalog, opy = 990)
}
\references{
  Hosking, J. R. M. and Wallis, J. R. (1995). A comparison of unbiased
  and plotting-position estimators of L moments. Water Resources
  Research. \bold{31}(8): 2019--2025.
  
  Ferro, C. and Segers, J. (2003). Inference for clusters of extreme
  values. Journal of the Royal Statistical Society B. \bold{65}:
  545--556.  
}
\author{Mathieu Ribatet}
\keyword{hplot}
