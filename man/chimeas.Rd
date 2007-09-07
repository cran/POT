\name{chimeas}
\alias{chimeas}
\title{Dependence Measures For Extreme Values Analysis}
\description{
  Provide two measures to assess for asymptotic dependence or independence
}
\usage{
chimeas(data, u.range, n.u = 500, xlab, ylabs, ci = 0.95,  boot = FALSE,
n.boot = 250, block.size = 50, show.bound = TRUE, which = 1:2, ask =
nb.fig < length(which) && dev.interactive(), \dots, col.ci = "grey",
col.bound = "blue", lty.ci = 1, lty.bound = 1)
}
\arguments{
  \item{data}{A matrix with 2 columns given the data.}
  \item{u.range}{Numeric vection of length 2 (may be missing): the range
    for the probabilities.}
  \item{n.u}{The number of probabilities to be considered}
  \item{xlab,ylabs}{The x-axis and ylabs labels. ylabs must be of length
    2}
  \item{ci}{The probability level for the confidence intervals}
  \item{boot}{Logical. If \code{TRUE}, confidence intervals are computed
    by bootstraping contiguous blocks. This may be needed if there is
    dependence between observations. If \code{FALSE} (the default),
    confidence intervals are derived using the Delta method.}
  \item{n.boot}{The number of bootstrap replicates.}
  \item{block.size}{The size of the ``contiguous'' blocks. See details.}
  \item{show.bound}{Logical. If \code{TRUE} (the default), the
    theoretical bound for the two statistics are plotted.}
  \item{which}{Which plot should be plotted? \code{1} for the
    \eqn{\chi}{chi} \code{2} for the \eqn{\overline{\chi}}{chibar}
    statistic and \code{1:2} for both of them.}
  \item{ask}{Logical. Should user be asked before each plot is
    computed?}
  \item{\dots}{Additional options to be passed to the \code{\link{plot}}
  function.}
  \item{col.ci,col.bound}{The color for the confidence
    intervals and theoretical bounds.}
  \item{lty.ci,lty.bound}{The line type for the confidence intervals and
    theoretical bounds.}
}
\value{
  A graphic window.
}
\details{
  These two plots help us to understand the dependence relationship
  between the two data set. The sign of \eqn{\chi(u)}{chi(u)} determines
  if the variables are positively or negatively correlated. Two variable
  are asymptotically independent if \eqn{\lim_{u\rightarrow1} \chi(u) =
    0}{lim_{u ~ 1} chi(u) = 0}. For the independent case, \eqn{\chi(u) =
    0}{chi(u) = 0} for all u in (0,1). For the perfect dependence case,
  \eqn{\chi(u) = 1}{chi(u) = 1} for all u in (0,1). Note that for a
  bivariate extreme value model, \eqn{\chi(u) = 2(1 - A(0.5))}{chi(u) =
    2(1-A(0.5))} for all u in (0,1).

  The measure \eqn{\overline{\chi}}{chibar} is only useful for
  asymptotically independent variables. Indeed, for asymptotically
  dependent variable, we have \eqn{\lim_{u\rightarrow
      1}\overline{\chi}(u) = 1}{lim_{u ~ 1}chibar(u) = 1}. For
  asymptotically independent variables, \eqn{\lim_{u\rightarrow
      1}\overline{\chi}(u)}{lim_{u ~ 1}chibar(u)} reflects the strength
  of the dependence between variables. For independent variables,
  \eqn{\overline{\chi}(u) = 0}{chibar(u)=0} for all u in (0,1).

  If there is (short range) dependence between observations, users may
  need to use bootstrap confidence intervals. Bootstrap series are
  obtained by sampling contiguous blocks, of length \code{l} say,
  uniformly with replacement from the original observations. The block
  length \code{l} should be chosen to be much greater than the
  short-range dependence and much smaller than the total number of
  observations.
}
\references{
  Coles, S., Heffernan, J. and Tawn, J. (1999) Dependence measures for
  extreme value analyses. \emph{Extremes} \bold{2} 339--365.
}
\seealso{
  \code{\link{tailind.test}}, \code{\link{specdens}},
  \code{\link{tsdep.plot}}
}
\examples{
mc <- simmc(2000, alpha = 0.9)
mc2 <- simmc(1000, alpha = 0.2)
##An independent case
par(mfrow = c(1,2))
chimeas(cbind(mc[1:1000], mc2))
##Asymptotic dependence
par(mfrow = c(1,2))
chimeas(cbind(mc[seq(1,2000, by = 2)], mc[seq(2,2000,by = 2)]))
##The same but with bootstrap ci
par(mfrow = c(1,2))
chimeas(cbind(mc[seq(1,2000, by = 2)], mc[seq(2,2000,by = 2)]), boot =
TRUE)
}
\author{Mathieu Ribatet}
\keyword{htest}


