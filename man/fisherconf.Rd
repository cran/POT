\name{Fisher Confidence Interval}
\alias{gpd.fishape}
\alias{gpd.fiscale}
\alias{gpd.firl}

\title{Fisher Based Confidence Interval for the GP Distribution}
\description{Compute Fisher based confidence intervals on parameter and
  return level for the GP distribution.  This is achieved through
  asymptotic theory and the Observed information matrix of Fisher and
  eventually the Delta method.
}

\usage{
gpd.fishape(object, conf = 0.95)
gpd.fiscale(object, conf = 0.95)
gpd.firl(object, prob, conf = 0.95)
}
\arguments{
  \item{object}{\code{R} object given by function \code{\link{fitgpd}}.}
  \item{prob}{The probability of non exceedance.}
  \item{conf}{Numeric. The confidence level.}
}
\value{
  Returns a vector of the lower and upper bound for the confidence
  interval.
}
\seealso{\code{\link{rp2prob}}, \code{\link{prob2rp}},
  \code{\link{gpd.pfscale}},
  \code{\link{gpd.pfshape}}, \code{\link{gpd.pfrl}} and
  \code{\link{confint}}
}
\author{Mathieu Ribatet}
\examples{
data(ardieres)
ardieres <- clust(ardieres, 4, 10 / 365, clust.max = TRUE)
f1 <- fitgpd(ardieres[,"obs"], 5, 'mle')
gpd.fishape(f1)
gpd.fiscale(f1)
}
\keyword{htest}
