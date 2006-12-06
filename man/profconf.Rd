\name{Confidence interval}
\alias{gpd.pfshape}
\alias{gpd.pfscale}
\alias{gpd.pfrl}

\title{Confidence interval for the GP distribution}
\description{Compute confidence intervals on parameter and return level
  for the GP distribution.  This is achieved through the profile
  likelihood procedure.
}

\usage{
gpd.pfshape(fitted, range, xlab, ylab, conf = 0.95, nrang = 100,
vert.lines = TRUE, \dots)
gpd.pfscape(fitted, range, xlab, ylab, conf = 0.95, nrang = 100,
vert.lines = TRUE, \dots)
gpd.pfrl(fitted, prob, range, thresh, xlab, ylab, conf = 0.95, nrang =
100, vert.lines = TRUE, \dots)
}
\arguments{
  \item{fitted}{\code{R} object given by function \code{\link{fitgpd}}.}
  \item{prob}{The probability of non exceedance.}
  \item{range}{Vector of dimension two. It gives the lower and upper
    bound on which the profile likelihood is performed.}
  \item{thresh}{Optional. The threshold. Only needed with non constant
    threshold.}
  \item{xlab, ylab}{Optional Strings. Allows to label the x-axis and
    y-axis. If missing, default value are considered.}
  \item{conf}{Numeric. The confidence level.}
  \item{nrang}{Numeric. It specifies the number of profile likelihood
    computed on the whole range \code{range}.}
  \item{vert.lines}{Logical. If \code{TRUE} (the default), vertical
    lines are plotted.}
  \item{\dots}{Optional parameters to be passed to the
    \code{\link{plot}} function.}  
}
\value{
  Returns a vector of the lower and upper bound for the profile confidence
  interval.  Moreover, a graphic of the profile likelihood function is
  displayed.
}

\references{
  Coles, S. (2001). \emph{An Introduction to Statistical
    Modelling of Extreme Values}. Springer Series in Statistics. London.
}
\examples{
data(ardieres)
fitted <- fitgpd(ardieres[, "flows"], 5, 'mle')
gpd.pfshape(fitted, c(0, 0.8))
rp2prob(10, 2)
gpd.pfrl(fitted, 0.95, c(14, 32))
}
\author{Mathieu Ribatet}
\keyword{htest}
