\name{Threshold Selection}
\alias{mrlplot}
\title{Threshold Selection: The Empirical Mean Residual Life Plot}
\description{
  The empirical mean residual life plot. 
}
\usage{
mrlplot(data, u.range, main, xlab, ylab, nt = max(100, length(data)),
lty = c(3, 1, 3), col = c('grey', 'black', 'grey'), conf = 0.95, lwd = c(1,
1.5, 1), \dots) 
}
\arguments{
  \item{data}{A numeric vector.}
  \item{u.range}{A numeric vector of length two, giving the limits for
    the thresholds at which the mean residual life plot is
    evaluated. If \code{u.range} is not given, sensible defaults
    are used.}
  \item{main}{Plot title.}
  \item{xlab, ylab}{x and y axis labels.}
  \item{nt}{The number of thresholds at which the mean residual life
    plot is evaluated.}
  \item{lty, col, lwd}{Arguments passed to \code{matplot}. The first
    and last elements of \code{lty} correspond to the lower and
    upper confidence limits respectively. Use zero to supress.}
  \item{conf}{The (pointwise) confidence coefficient for the
    plotted confidence intervals.} 
  \item{\dots}{Other arguments to be passed to \code{matplot}.}
}
\details{
  The empirical mean residual life plot is the locus of points
  \deqn{\left(u,\frac{1}{n_u} \sum\nolimits_{i=1}^{n_u}
    (x_{(i)} - u) \right)}{{u,1/n_u \sum_{i=1}^{n_u} (x(i) - u)}}
  where \eqn{x_{(1)}, \dots, x_{(n_u)}}{x(1), \dots, x(n_u)} are
  the \eqn{n_u} observations that exceed the threshold \eqn{u}.  If the
  exceedances of a threshold \eqn{u_0}{u0} are generalized Pareto, the
  empirical mean residual life plot should be approximately linear for
  \eqn{u > u_0}{u > u0}.
  
  The confidence intervals within the plot are symmetric intervals
  based on the approximate normality of sample means.
}
\value{
  A list with components \code{x} and \code{y} is invisibly returned.
  The components contain those objects that were passed to the formal
  arguments \code{x} and \code{y} of \code{matplot} in order to create
  the mean residual life plot.
}
\author{Stuart Coles and Alec Stephenson}  
\seealso{\code{\link{fitgpd}}, \code{\link{matplot}},
  \code{\link{tcplot}}}

\references{Coles, S. (2001) \emph{An Introduction to Statistical
    Modelling of Extreme Values}. Springer Series in Statistics. London.

  Embrechts, P., Kl\"uppelberg, C., and Mikosch, T. (1997)
  \emph{Modelling Extremal Events for Insurance and Finance.}
}
\examples{
data(ardieres)
ardieres <- clust(ardieres, 4, 10 / 365, clust.max = TRUE)
flows <- ardieres[, "obs"]
mrlplot(flows)
}
\keyword{hplot}

