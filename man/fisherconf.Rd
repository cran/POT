\name{Confidence interval}
\alias{gpd.fishape}
\alias{gpd.fiscale}
\alias{gpd.firl}

\title{Confidence interval for the GP distribution}
\description{Compute confidence intervals on parameter and return level
  for the GP distribution. This is achieved through asymptotic theory
  and the Observed information matrix of Fisher.
}

\usage{
gpd.fishape(fitted, conf = 0.95)
gpd.fiscale(fitted, conf = 0.95)
gpd.firl(fitted, prob, conf = 0.95)
}
\arguments{
  \item{fitted}{\code{R} object given by function \code{\link{fitgpd}}.}
  \item{prob}{The probability of non exceedance.}
  \item{conf}{Numeric. The confidence level.}
}
\value{
  Returns a vector of the lower and upper bound for the confidence
  interval.
}
\seealso{\code{\link{rp2prob}} and \code{\link{prob2rp}}.}
\author{Mathieu Ribatet}
\examples{
data(ardieres)
fitted <- fitgpd(ardieres[,"flows"], 5, 'mle')
gpd.fishape(fitted)
gpd.fiscale(fitted)
}
\keyword{htest}
