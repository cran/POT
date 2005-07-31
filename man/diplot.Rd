\name{Threshold Selection}
\alias{diplot}

\title{Dispersion Index Plot}
\description{The Dispersion Index Plot}
\usage{
diplot(data, tlim, conf=0.95, main = 'Dispersion Index Plot',\dots)
}
\arguments{
  \item{data}{A matrix with two column. The first one represents the
    date of events and the second the data associated with those dates.}
  \item{tlim}{A numeric vector of length two giving the limit of
    threshold analyzed. If missing, default values are taken.}
  \item{conf}{The confident coefficient for the plotted confidence intervals.}
  \item{main}{The title of the plot.}
  \item{...}{Other arguments to be passed to the \code{plot} function.}
}
\value{It returns invisibly a list with two components. The first one
  \code{'thresh'} gives the thresholds analyzed. The second \code{'DI'}
  gives the dispersion index relative to the threshold.
}
\details{
  According to the Extreme Value Theory, the number of exceedance
  over a high threshold in a fixed period - generally a year - must be
  distributed as Poisson processus. As for a random variable Poisson
  distributed, the ratio of the variance and the mean is equal to 1, one
  can test if the ratio \eqn{\code{DI} = var / mean} differs from
  1. Moreover, confidence levels for \code{DI} can be calculated by
  testing against a \eqn{\chi^2} distribution with \code{M}-1 degree of
  freedom, \code{M} being the total number of fixed periods -generally
  the total number of years in the sample. So, the Poisson hypothesis is
  not rejected if the estimated \code{DI} is within the range
  \deqn{\left[ \frac{\chi^2_{\alpha/2, \code{M}-1}}{\code{M}-1},
    \frac{\chi^2_{1 - \alpha/2, \code{M}-1} }{\code{M} - 1} \right]}{%
    \left[\chi^2_{\alpha/2, \code{M}-1} / ( \code{M} - 1 ), \chi^2_{1
      -\alpha/2, \code{M}-1} / ( \code{M} - 1 ) \right]}
}
\examples{
data(ardieres)
diplot(ardieres)
}
\keyword{hplot}