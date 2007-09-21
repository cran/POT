\name{fitexi}
\alias{fitexi}
\title{Extremal Index Estimation}
\description{
  Estimation of the extremal index using interexceedance times.
}
\usage{
fitexi(data, threshold)
}
\arguments{
  \item{data}{A matrix with two columns: \code{``obs''} for the
    observations and \code{``time''} for the time.}
  \item{threshold}{The threshold.}
}
\value{
  This function returns a list with two components. The first one
  \code{``exi''} gives the estimate of the extremal index; while the
  second, \code{``tim.cond''} gives the time condition for independence
  between events to be passed to the \code{\link{clust}} function.
}
\details{
  The extremal index estimator proposed by Ferro and Segers
  (2003) is based on interexceedance times. In particular, it does not
  require a specific declusterization of the time series.

  The \code{``tim.cond''} gives an ``automatic'' procedure to decluster
  the time series without any subjective choice to define the
  independence condition between clusters.
}
\references{
  Ferro, C. and Segers, J. (2003) Inference for clusters of extreme
  values. Journal of the Royal Statistical Society. Series B \bold{65:2}
  545--556.
}
\seealso{
  \code{\link{clust}}
}
\examples{
data(ardieres)
exi <- fitexi(ardieres, 6)
exi
##Extract events with the ``automatic'' procedure
##Please be carefull for the time unit in tim.cond
clust(ardieres, 6, tim.cond = exi$tim.cond / 365.25, clust.max = TRUE)
}
\author{Mathieu Ribatet}
\keyword{htest}
