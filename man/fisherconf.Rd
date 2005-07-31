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
gpd.firl(fitted, retper, mu, conf = 0.95)
}
\arguments{
  \item{fitted}{\code{R} object given by function \code{\link{fitgpd}}.}
  \item{retper}{. The return period.}
  \item{mu}{Numeric. The mean number of event per year.}
  \item{conf}{Numeric. The confidence level.}
}
\value{
  Returns a vector of the lower and upper bound for the confidence
  interval.
}
\examples{
data(ardieres)
fitted <- fitgpd(ardieres[,"flows"], 5, 'mle')
gpd.fishape(fitted)
gpd.fiscale(fitted)
}
\keyword{htest}