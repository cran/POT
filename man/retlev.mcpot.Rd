\name{retlev.mcpot}
\alias{retlev.mcpot}
\title{Return Level Plot: Markov Chains Models}
\description{
  Return level plot for Markov chain POT models. 
}
\usage{
\method{retlev}{mcpot}(fitted, opy, exi, main, xlab, ylab, xlimsup,
\dots)
}
\arguments{
  \item{fitted}{A object of class \code{``uvpot''}. Most often, the
    return of the \code{\link{fitgpd}} function.}
  \item{opy}{The number of \bold{O}bservations \bold{P}er \bold{Y}ear (or
    more generally per block). If missing, it is set it to 365 i.e. daily
    values with a warning.}
  \item{exi}{Numeric. The extremal index. If missing, an estimate is
    given using the \code{\link{fitexi}} function.}
  \item{main}{The title of the graphic. If missing, the title is set to
    ``Return Level Plot''.}
  \item{xlab,ylab}{The labels for the x and y axis. If missing, they are
    set to ``Return Period (Years)'' and ``Return Level'' respectively.}
  \item{xlimsup}{Numeric. The right limit for the x-axis. If missing, it
    is setted to 500.}
  \item{\dots}{Other arguments to be passed to the \code{\link{plot}}
    function.}
}
\value{
  A graphical window. In addition, it returns invisibly the return level
  function.
}
\details{
  Let \eqn{X_1, \ldots,X_n}{X_1, \ldots,X_n} be the first \eqn{n}{n}
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
  Though this is computationally expensive, we recommend to give the
  extremal index estimate using the \code{\link{dexi}} function. Indeed,
  there is a severe bias when using the Ferro and Segers (2003)
  estimator - as it is estimated using observation and not the Markov
  chain model.
}
\examples{
data(ardieres)
Mcalog <- fitmcgpd(ardieres[,"obs"], 5, "alog")
retlev(Mcalog, opy = 990)
}
\seealso{
  \code{\link{retlev}}, \code{\link{retlev.uvpot}},
  \code{\link{retlev.mcpot}}, \code{\link{fitexi}}
}
\references{
  Ferro, C. and Segers, J. (2003). Inference for clusters of extreme
  values. Journal of the Royal Statistical Society B. \bold{65}:
  545--556.
}
\author{Mathieu Ribatet}
\keyword{hplot}
